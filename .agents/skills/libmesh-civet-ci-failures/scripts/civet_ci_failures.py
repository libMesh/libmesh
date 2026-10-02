#!/usr/bin/env python3
"""
Report CIVET failures for a libMesh pull request, commit or job in a compact, bounded form.

Reads only through the GitHub API (via the gh CLI) and CIVET's public job logs,
so it needs no CIVET credentials. Two sources are combined:

  - The commit statuses CIVET posts per job, which give job-level state and
    the URL of each CIVET job.
  - The step logs of each failed job, which give the build diagnostics, the
    failing unit tests and examples, and how to reproduce them.

There are three ways to name what to report on. --pr covers a pull request,
--sha covers a commit, which is what makes pushes to devel reportable, and
--job covers a single CIVET job named by its URL. --job reports the same
findings as the other two and additionally resolves the URL to the event,
recipe and commits behind it, none of which the URL itself carries.

A libMesh check step runs `make check`: the contrib test suites, the CppUnit
unit tests, and every example. A step CIVET killed at its time limit leaves no
failure message, so for those the unit test or example that was still running
is reported instead, with the unit test identified by its position in the run.

The output is deliberately small. Passing tests are never reported and every
list is capped. Use --json for the merged payload.

Adapted from MOOSE's python/civet_ci_failures/civet_ci_failures.py, which
reads the same CIVET data for the MOOSE TestHarness.
"""

import argparse
import collections
import configparser
import glob
import json
import os
import re
import subprocess
import sys
from typing import List, NamedTuple, Optional, Sequence

# Remotes tried, in order, when none is named. A developer checkout usually
# has origin pointing at a personal fork rather than the upstream, so origin
# is the last resort. Override with --upstream-remote, --repo, or the
# CIVET_UPSTREAM_REMOTE environment variable.
UPSTREAM_REMOTE_CANDIDATES = ("up", "upstream", "origin")

# Default cap on the number of failing tests to describe
DEFAULT_MAX_FAILURES = 20

# Default cap on the number of failed jobs to list
MAX_JOBS = 15

# Status state CIVET posts for a job that has not finished. GitHub's combined
# state cannot be used to tell whether an event is done: it turns "failure" as
# soon as any one job fails, with every other job still pending.
PENDING_STATE = "pending"

# Description CIVET sets on a job that was skipped because something it
# depends on failed. These are cascade effects, not failures to fix.
BLOCKED_DESCRIPTION = "Won't run due to failed dependencies"

# Every CIVET step log ends with this, which is how a failing step is found
# without asking CIVET for per-step status
STEP_RETURN_CODE_RE = re.compile(r"completed with return code (\d+)")

# A compiler or linker diagnostic. The make failure lines that follow are
# noise once the diagnostic itself is in hand, so they are only a fallback.
DIAGNOSTIC_RE = re.compile(r"\berror:|\bundefined reference to\b|^ld: ")
MAKE_FAILURE_RE = re.compile(r"\*\*\* \[.*\] Error \d+")

# Failures from the surrounding machinery rather than from compiled code or a
# test: a container build, a download, or a script the recipe drives. These
# carry no compiler diagnostic, so they need their own patterns.
INFRA_ERROR_RE = re.compile(
    r"^FATAL:|curl: \(\d+\)|CalledProcessError|Connection reset by peer|"
    r"Could not resolve host|No space left on device|Connection timed out|"
    r"^fatal: unable to access|^fatal: remote error|Empty reply from server|"
    r"RPC failed"
)

# Lines that match the diagnostic pattern while only reporting someone else's
# failure. make relaying an error from a recipe it ran, or git reporting a
# transport error, says nothing about the cause, so these must not preempt the
# infra check that finds it.
NON_DIAGNOSTIC_RE = re.compile(r"^(?:make|gmake)(?:\[\d+\])?: |^fatal: ")

# Step teardown, which is the last thing in every log and says nothing about
# why the step failed. It has to be dropped before falling back to the tail.
TEARDOWN_MARKERS = (
    "Removing ",
    "rm -rf",
    "Submitting step statistics",
    "Failed to submit step statistics",
    "Execution in ",
    "Collecting test statistics",
    "ERROR: Exiting with code",
)

# The container a step executed in, named once in that step's header
CONTAINER_RE = re.compile(r"Executing \S+ in (\S+://\S+|\S*\.sif)")

# An assignment in the CIVET environment that every step dumps in its header.
# This is the only place a job URL alone can be resolved to the event, commit
# and recipe it belongs to, none of which are in the URL.
CIVET_ENV_RE = re.compile(r'^(CIVET_[A-Z_]+)="(.*)"$', re.MULTILINE)

# CIVET's notice that it killed a step at the recipe's time limit. Everything
# printed after it is the kill propagating (PETSc's signal handler reporting
# SIGTERM, MPI_Abort, make reporting Terminated), not the cause.
CANCEL_RE = re.compile(
    r"CIVET: Cancelling job due to step taking longer than the max (\d+) seconds"
)

# The unit test invocation `make check` logs, with or without an MPI launcher
# (LIBMESH_RUN) in front of it
UNIT_RUN_RE = re.compile(
    r"^(?P<launcher>(?:\S*/)?(?:mpiexec|mpirun|srun)\b.*?\s)?"
    r"\./(?P<binary>unit_tests-\w+)(?P<args>\s.*)?$"
)

# CppUnit's end-of-run results: the success line, the failure tally, and the
# header of each failure's detail, "1) test: Suite::test (F) line: 12 file.C",
# where F marks a failed assertion and E an uncaught exception
CPPUNIT_OK_RE = re.compile(r"^OK \((\d+) tests?\)")
CPPUNIT_TALLY_RE = re.compile(r"^Run:\s+(\d+)\s+Failures:\s+(\d+)\s+Errors:\s+(\d+)")
CPPUNIT_FAILURE_RE = re.compile(
    r"^(\d+)\) test: (\S+) \(([FE])\)(?: line: (\d+) (\S+))?"
)

# Lines of a CppUnit failure detail that only restate that a failure happened
CPPUNIT_DETAIL_BOILERPLATE = ("assertion failed", "equality assertion failed")

# The banners examples/run_common.sh prints around each example run, each
# followed by a line carrying the command. A banner with no closing one is an
# example that exited nonzero, or was still running when the step was killed.
EXAMPLE_START_RE = re.compile(r"^\* Running Example (\S+):$")
EXAMPLE_DONE_RE = re.compile(r"^\* Done Running Example (\S+):$")
EXAMPLE_COMMAND_RE = re.compile(r"^\*\s+(\S.*)$")

# The directory make enters for an example, which is where its command runs
EXAMPLE_DIR_RE = re.compile(r"Entering directory '[^']*/(examples/[^']+)'")

# Lines that state why a libMesh program died: a failed libmesh_assert, a
# libmesh_error or uncaught exception, a fatal signal, or an error a script
# reported. PETSc's own boilerplate around a signal is not matched.
ERROR_LINE_RE = re.compile(
    r"Assertion `.*' failed|^terminate called|what\(\): |libMesh terminating|"
    r"Caught signal number \d+|^(?:\[\d+\])?ERROR\b|^[A-Za-z_:]*Error: "
)

# Applied before comparing two error messages, so that a differing temporary
# directory or index does not make one failure look like many
NORMALIZE_SUBS = (
    (re.compile(r"/tmp/[^/\s]+"), "/tmp/<dir>"),
    (re.compile(r"0x[0-9a-fA-F]+"), "<addr>"),
    (re.compile(r"\d+"), "N"),
)

# Status labels for the failures found in a check step
HANG = "HANG"
CRASH = "CRASH"
FAILURE = "FAILURE"
ERROR = "ERROR"

# Known remediations, tried in order against a failure label or error
# signature
REMEDIATION_HINTS = (
    (
        HANG,
        "nothing failed visibly; CIVET killed the step at its time limit. If the "
        "hang appears only at higher process counts, suspect an assertion on "
        "rank-local data that fails on some ranks, which leaves the others "
        "blocked in the next collective call",
    ),
    (
        CRASH,
        "the program ended without a CppUnit result; look for a signal or an "
        "uncaught error earlier in the same run",
    ),
    (
        "curl: (",
        "a download failed during the build; spurious, re-run the job rather "
        "than changing code",
    ),
    (
        "fatal: unable to access",
        "a git clone or fetch failed during the build, so the remote was never "
        "reached; spurious as far as the branch is concerned",
    ),
    (
        "Could not resolve host",
        "DNS failed for the host being cloned from; spurious, re-run once the "
        "host is reachable",
    ),
    (
        "fatal: remote error",
        "the ref being fetched is not on the remote, usually a submodule bumped "
        "to a commit that has not been pushed upstream; push it, then re-run",
    ),
    (
        "FATAL:",
        "the container the step runs in could not be pulled or started; "
        "suspect infrastructure before code, and re-run",
    ),
)

# Default cap on error signatures reported per job
DEFAULT_MAX_SIGNATURES = 4

# Default cap on unique diagnostics reported per failing step
DEFAULT_MAX_DIAGNOSTICS = 5

# Default cap on jobs to pull logs for, since each is a tarball download
DEFAULT_MAX_LOG_JOBS = 6

# Default cap on lines printed from a single step's log. A step can be several
# megabytes, so raw log access is always bounded.
DEFAULT_MAX_LOG_LINES = 100


class GitHubError(SystemExit):
    """Exception for a failed gh invocation."""


def gh(args: Sequence[str]) -> str:
    """Run a gh command and return its stdout."""
    try:
        result = subprocess.run(
            ["gh", *args], capture_output=True, text=True, check=False
        )
    except FileNotFoundError as e:
        raise GitHubError("The gh CLI is required and was not found on PATH") from e
    if result.returncode != 0:
        raise GitHubError(f"gh {' '.join(args)} failed:\n{result.stderr.strip()}")
    return result.stdout


def gh_json(args: Sequence[str]):
    """Run a gh command that produces JSON and parse it."""
    out = gh(args).strip()
    return json.loads(out) if out else None


def latest_statuses(statuses: Sequence[dict]) -> List[dict]:
    """
    Reduce the status list to the newest status per context.

    GitHub returns every status ever posted for a commit, newest first, and
    CIVET posts one when a job starts and again when it finishes. Taking the
    first occurrence of each context avoids reporting a finished job as
    pending.
    """
    seen = set()
    latest = []
    for status in statuses:
        context = status.get("context")
        if context in seen:
            continue
        seen.add(context)
        latest.append(status)
    return latest


def current_branch() -> str:
    """Get the name of the checked out branch."""
    result = subprocess.run(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise GitHubError("Not in a git repository; pass --pr")
    return result.stdout.strip()


def repo_from_remote(remote: str) -> Optional[str]:
    """Get the owner/name slug named by a git remote, if it exists."""
    result = subprocess.run(
        ["git", "remote", "get-url", remote],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        return None
    match = re.search(r"[:/]([^/:]+)/([^/]+?)(?:\.git)?$", result.stdout.strip())
    return f"{match.group(1)}/{match.group(2)}" if match else None


def resolve_repo(repo: Optional[str], remote: Optional[str]) -> Optional[str]:
    """
    Decide which repository to query.

    Precedence: an explicit --repo, then a named remote (--upstream-remote or
    CIVET_UPSTREAM_REMOTE), then each candidate remote in turn.
    """
    if repo:
        return repo
    named = remote or os.environ.get("CIVET_UPSTREAM_REMOTE")
    for candidate in [named] if named else UPSTREAM_REMOTE_CANDIDATES:
        if slug := repo_from_remote(candidate):
            return slug
    return None


def resolve_target(repo: str, pr: Optional[int], sha: Optional[str]) -> dict:
    """
    Resolve what to report on: a pull request, or a bare commit.

    CIVET attaches its statuses to the commit rather than to the pull request,
    so a SHA is sufficient on its own. That is what lets this report on a push
    to devel, where there is no pull request at all.
    """
    if sha:
        return {"number": None, "headRefName": None, "headRefOid": sha, "url": None}
    return resolve_pr(repo, pr)


def resolve_pr(repo: str, pr: Optional[int]) -> dict:
    """
    Resolve the PR number and head SHA, defaulting to the branch's PR.

    The repository is always passed explicitly. A libMesh checkout commonly has
    a remote per fork whose PRs are reviewed locally, and gh picks one of them
    arbitrarily when it has no default, which resolves to the wrong repository
    or to no pull request at all.
    """
    fields = "number,headRefOid,headRefName,url"

    if pr is not None:
        info = gh_json(["pr", "view", str(pr), "--repo", repo, "--json", fields])
        if not info:
            raise GitHubError(f"Could not resolve PR #{pr} in {repo}")
        return info

    branch = current_branch()
    found = gh_json(
        [
            "pr",
            "list",
            "--repo",
            repo,
            "--head",
            branch,
            "--json",
            fields,
            "--limit",
            "1",
        ]
    )
    if not found:
        raise GitHubError(
            f"No open pull request in {repo} for branch '{branch}'; pass --pr"
        )
    return found[0]


def fetch_state(slug: str, sha: str) -> dict:
    """
    Fetch the rollup state and every per-job status for a SHA.

    The combined status endpoint caps its statuses array at 30, which silently
    truncates an event with more jobs than that, so the per-job statuses come
    from the paginated list endpoint instead. Its rollup state is still taken
    from the combined endpoint, where it accounts for every job.
    """
    combined = gh_json(["api", f"repos/{slug}/commits/{sha}/status"]) or {}

    # --jq emits one object per line, avoiding the concatenated arrays that
    # --paginate produces on its own
    raw = gh(
        [
            "api",
            f"repos/{slug}/commits/{sha}/statuses",
            "--paginate",
            "--jq",
            ".[]",
        ]
    )
    statuses = [json.loads(line) for line in raw.splitlines() if line.strip()]

    return {
        "state": combined.get("state", "unknown"),
        "statuses": latest_statuses(statuses),
    }


def fetch_job_steps(job_url: str) -> List[tuple]:
    """
    Fetch the step logs for a CIVET job as (name, text) pairs.

    CIVET serves every step of a job as a gzipped tarball with one member per
    step, which needs no authentication for public repositories. The tarball is
    read straight into memory and nothing is written to disk, so there is no
    temporary state to clean up; the cost is that a job's whole log is resident
    while it is being parsed.
    """
    import io
    import tarfile
    import urllib.request

    base = job_url.rstrip("/").replace("/job/", "/job_results/")
    with urllib.request.urlopen(f"{base}/") as response:
        payload = response.read()

    steps = []
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as tar:
        for member in tar.getmembers():
            handle = tar.extractfile(member)
            if handle is None:
                continue
            text = handle.read().decode("utf-8", "replace")
            steps.append((member.name.split("/")[-1], text))
    return steps


def extract_job_info(steps: Sequence[tuple]) -> dict:
    """
    Read the CIVET environment that a job's steps dump in their headers.

    A job URL says nothing on its own about what the job was testing. The
    recipe name, the pull request and both commits live only in this dump, and
    they are what decide whether a failure belongs to the branch under test or
    to the base. Later steps overwrite earlier ones, which only matters for the
    per-step entries; everything used here is a property of the job.
    """
    env: dict = {}
    for _, text in steps:
        env.update(CIVET_ENV_RE.findall(text))
    return env


def print_job_header(job_url: str, env: dict) -> None:
    """Print the event, recipe and commits behind a single CIVET job."""
    print(f"{env.get('CIVET_RECIPE_NAME') or '?'}  {job_url}")

    if pr := env.get("CIVET_PR_NUM"):
        print(f"  pull request {pr} ({env.get('CIVET_EVENT_CAUSE') or '?'})")
    elif cause := env.get("CIVET_EVENT_CAUSE"):
        print(f"  {cause}")

    head_repo = env.get("CIVET_HEAD_REPO") or "?"
    head_ref = env.get("CIVET_HEAD_REF") or "?"
    head_sha = (env.get("CIVET_HEAD_SHA") or "?")[:12]
    print(f"  head {head_repo}:{head_ref} @ {head_sha}")

    base_repo = env.get("CIVET_BASE_REPO") or "?"
    base_ref = env.get("CIVET_BASE_REF") or "?"
    base_sha = (env.get("CIVET_BASE_SHA") or "?")[:12]
    print(f"  base {base_repo}:{base_ref} @ {base_sha}")

    # An invalidated job reran on the same commit, so its logs may predate the
    # state the rest of the event was built from
    if env.get("CIVET_INVALIDATED") == "True":
        print("  this job was invalidated and rerun")


class Failure(NamedTuple):
    """One occurrence of a failing test: where it failed and how to repeat it."""

    # Display name of the job it failed in
    job: str
    # Command that reruns just this test in the mode that failed it
    command: str
    # The mode itself, ignoring resource limits, for counting distinct modes
    mode: str
    # Container the step ran in, if its header named one
    container: Optional[str]


def unique(values: Sequence[str]) -> List[str]:
    """Deduplicate while preserving order."""
    seen = set()
    out = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def extract_container(text: str) -> Optional[str]:
    """
    Get the container the given step executed in, if its header names one.

    Takes a single step's log, because the container is a property of the step
    and not of the job: a fetch step runs in a base image while the steps that
    build and test run in a versioned one. That is why the caller pairs what
    this returns with a particular step's invocation.
    """
    matches = CONTAINER_RE.findall(text)
    return matches[-1] if matches else None



def parse_runs(lines: Sequence[str]) -> List[dict]:
    """
    Split a check step's log into its unit test and example runs.

    Each run records the span of the log it owns and whether it reached its own
    conclusion: a CppUnit result for the unit tests, the closing banner for an
    example. A run that did not either exited abnormally or was still running
    when the step was killed. A unit test run keeps its span past its result,
    up to the next run, because CppUnit prints the failure details after it.
    """
    runs: List[dict] = []
    current: Optional[dict] = None
    example_dir = None

    def close(end: int) -> None:
        nonlocal current
        if current is not None:
            current["end"] = end
            runs.append(current)
            current = None

    for i, line in enumerate(lines):
        stripped = line.strip()

        if match := EXAMPLE_DIR_RE.search(stripped):
            example_dir = match.group(1)
        elif match := UNIT_RUN_RE.match(stripped):
            close(i)
            launcher = (match.group("launcher") or "").strip()
            args = (match.group("args") or "").strip()
            current = {
                "kind": "unit",
                "name": match.group("binary"),
                "start": i,
                "directory": "tests",
                "command": " ".join(
                    part for part in (launcher, f"./{match.group('binary')}", args) if part
                ),
                "finished": False,
            }
        elif match := EXAMPLE_START_RE.match(stripped):
            close(i)
            following = lines[i + 1].strip() if i + 1 < len(lines) else ""
            command = EXAMPLE_COMMAND_RE.match(following)
            current = {
                "kind": "example",
                "name": match.group(1),
                "start": i,
                "directory": example_dir,
                "command": command.group(1) if command else None,
                "finished": False,
            }
        elif current is None:
            continue
        elif current["kind"] == "unit" and (
            CPPUNIT_OK_RE.match(stripped) or CPPUNIT_TALLY_RE.match(stripped)
        ):
            current["finished"] = True
        elif current["kind"] == "example" and EXAMPLE_DONE_RE.match(stripped):
            current["finished"] = True
            close(i + 1)

    close(len(lines))
    return runs


def progress_marks(lines: Sequence[str]) -> tuple:
    """
    Count the tests CppUnit started, and how many of them failed, from its progress marks.

    CppUnit prints a dot as each test starts and an F or E as one fails, with no
    test names, so the count of dots is what identifies the test a run died in:
    with k dots, it is the kth test of the run. Output a test prints interleaves
    with the marks, so only marks at the start of a line are counted, and a
    trailing F or E that begins a word, as in "Failures", is output rather than
    a mark.
    """
    started = failed = 0
    for line in lines:
        match = re.match(r"[.FE]+", line)
        if not match:
            continue
        marks = match.group(0)
        if line[match.end() : match.end() + 1].isalnum():
            marks = marks.rstrip("FE")
        started += marks.count(".")
        failed += len(marks) - marks.count(".")
    return started, failed


def cppunit_failures(lines: Sequence[str]) -> List[dict]:
    """
    Get the failure details CppUnit printed after a run's tally.

    Each detail opens with a header naming the test, followed by lines saying
    what failed, up to a blank line. The lines that only say an assertion
    failed are dropped, so the expression or exception message is what remains.
    """
    failures = []
    for i, line in enumerate(lines):
        match = CPPUNIT_FAILURE_RE.match(line.strip())
        if not match:
            continue
        detail = []
        for following in lines[i + 1 :]:
            following = following.strip()
            if not following or CPPUNIT_FAILURE_RE.match(following):
                break
            if following not in CPPUNIT_DETAIL_BOILERPLATE:
                detail.append(following)
        location = f"{match.group(5)}:{match.group(4)}" if match.group(4) else None
        failures.append(
            {
                "test": match.group(2),
                "status": FAILURE if match.group(3) == "F" else ERROR,
                "detail": " / ".join(detail[:3]),
                "location": location,
            }
        )
    return failures


def normalize_error(message: str) -> str:
    """Reduce an error message to a form two occurrences can be compared on."""
    for pattern, replacement in NORMALIZE_SUBS:
        message = pattern.sub(replacement, message)
    return message


def error_signatures(lines: Sequence[str]) -> "collections.Counter":
    """
    Count the distinct root errors in some output.

    A single defect usually surfaces in many tests, so the count of unique
    signatures says how many things actually went wrong, which the count of
    failing tests does not.
    """
    signatures: collections.Counter = collections.Counter()
    for line in lines:
        if ERROR_LINE_RE.search(line.strip()):
            signatures[normalize_error(line.strip())] += 1
    return signatures


def check_failures(lines: Sequence[str], runs: Sequence[dict], hung: bool) -> List[dict]:
    """
    Describe each failing unit test and example of a check step.

    A unit test run reports its own failures. A run that never reached its
    conclusion is reported as a crash, or, when the step was killed and it is
    the last run, as a hang; for the unit tests the test is then named by its
    position, which is all the log records.
    """
    failures = []
    for index, run in enumerate(runs):
        span = lines[run["start"] : run["end"]]
        prefix = f"cd {run['directory']} && " if run["directory"] else ""
        command = f"{prefix}{run['command'] or '?'}"
        status = HANG if hung and index == len(runs) - 1 else CRASH

        if run["kind"] == "unit":
            for failure in cppunit_failures(span):
                failures.append(
                    {
                        "test": failure["test"],
                        "status": failure["status"],
                        "reason": failure["detail"]
                        + (f" ({failure['location']})" if failure["location"] else ""),
                        "command": f"{command} --re '^{failure['test']}$'",
                        "mode": run["command"],
                    }
                )
            if run["finished"]:
                continue
            started, failed = progress_marks(span)
            failures.append(
                {
                    "test": f"{run['name']} test #{started}",
                    "status": status,
                    "reason": f"{started} tests started, {failed} failed",
                    # Only a local run with test names printed can name it
                    "command": f"{command} --verbose",
                    "mode": run["command"],
                }
            )
        elif not run["finished"]:
            signatures = error_signatures(span)
            failures.append(
                {
                    "test": f"example {run['name']}",
                    "status": status if status == HANG else FAILURE,
                    "reason": signatures.most_common(1)[0][0] if signatures else None,
                    "command": command,
                    "mode": run["command"],
                }
            )
    return failures


def failing_steps(steps: Sequence[tuple]) -> List[tuple]:
    """
    Select the steps that exited nonzero, as (name, code, text).

    A job does not necessarily stop at its first failing step, and steps that
    run afterwards can succeed, so the last step is not a reliable indicator.
    """
    failing = []
    for name, text in steps:
        codes = STEP_RETURN_CODE_RE.findall(text)
        code = int(codes[-1]) if codes else None
        if code:
            failing.append((name, code, text))
    return failing


def extract_step_errors(text: str, limit: int) -> dict:
    """
    Extract the actionable lines from one failing step's log.

    A check step is read run by run, for failing unit tests and examples and
    for the one still running when CIVET killed the step. Anything else, or a
    check step none of whose runs failed, is read as a build: a compiler
    diagnostic is repeated once per translation unit and must be deduplicated.
    """
    lines = [line.rstrip() for line in text.splitlines()]

    cancel_at = next(
        (i for i, line in enumerate(lines) if CANCEL_RE.search(line)), None
    )
    limit_seconds = (
        int(CANCEL_RE.search(lines[cancel_at]).group(1)) if cancel_at is not None else None
    )

    # CppUnit's progress marks are flushed as the kill arrives, so the runs
    # span the whole log, but their error output stops at the kill
    runs = parse_runs(lines)
    failures = check_failures(lines, runs, cancel_at is not None)
    if failures:
        before_kill = lines[:cancel_at] if cancel_at is not None else lines
        signatures: collections.Counter = collections.Counter()
        for run in runs:
            if not run["finished"]:
                signatures.update(error_signatures(before_kill[run["start"] : run["end"]]))
        return {
            "kind": "check",
            "failures": failures,
            "killed_after": limit_seconds,
            "signatures": signatures,
        }

    def entry(kind: str, items: List[str]) -> dict:
        return {
            "kind": kind,
            "items": items[:limit],
            "all_items": items,
            "total": len(items),
            "killed_after": limit_seconds,
        }

    if cancel_at is not None:
        lines = lines[:cancel_at]

    # A real compiler or linker diagnostic names the cause and wins outright
    diagnostics = unique(
        [
            line
            for line in lines
            if DIAGNOSTIC_RE.search(line) and not NON_DIAGNOSTIC_RE.match(line.strip())
        ]
    )
    if diagnostics:
        return entry("build", diagnostics)

    # Otherwise the surrounding machinery is where the cause is, and only then
    # the build system's own complaint about a failure it inherited
    infra = unique([line for line in lines if INFRA_ERROR_RE.search(line.strip())])
    if infra:
        return entry("infra", infra)

    symptoms = unique([line for line in lines if MAKE_FAILURE_RE.search(line)])
    if symptoms:
        return entry("build", symptoms)

    # Nothing recognized, so fall back to the tail, less the teardown that
    # every log ends with and that never explains the failure
    tail = [
        line
        for line in lines
        if line.strip() and not any(m in line for m in TEARDOWN_MARKERS)
    ][-limit:]
    return entry("hang" if cancel_at is not None else "unknown", tail)


def hints_for(texts: Sequence[str]) -> List[str]:
    """Get the remediations that apply to the given labels and signatures."""
    hints = []
    for key, hint in REMEDIATION_HINTS:
        if any(key in text for text in texts):
            hints.append(f"{key} -> {hint}")
    return hints


def rollup_entries(rollup: dict, limit: int) -> List[dict]:
    """
    Summarize the failing tests, one entry per test.

    The reproduce command is the shortest of the modes that failed, since a
    longer one only adds variables that were not needed to fail. The mode
    count is what says whether a failure is specific to one way of running,
    such as one process count.
    """
    entries = []
    for test, info in list(rollup.items())[:limit]:
        failures = info["failures"]
        jobs = unique([f.job for f in failures])
        commands = unique([f.command for f in failures])
        modes = unique([f.mode for f in failures])
        # The container has to be the one belonging to the invocation reported
        # below, not just any the test failed under
        chosen = min(commands, key=len)
        container = next(
            (f.container for f in failures if f.command == chosen and f.container), None
        )
        reasons = sorted(info.get("reasons", ()))
        entries.append(
            {
                "test": test,
                "status": info["status"],
                "reasons": reasons,
                "hints": hints_for([info["status"], *reasons]),
                "jobs": jobs,
                "modes": len(modes),
                "reproduce": chosen,
                "container": container,
            }
        )
    return entries


def collect_job_errors(
    context: str,
    url: str,
    steps: Sequence[tuple],
    rollup: dict,
    hint_texts: set,
    max_diagnostics: int,
    max_signatures: int,
) -> dict:
    """
    Gather what one job's step logs contain, accumulating shared findings.

    The test rollup and the remediation hints span jobs, so they are passed in
    and added to rather than returned. Splitting this out is what lets a single
    job named by its URL and a whole event's failed jobs produce the same
    findings from the same code.
    """
    entry: dict = {"context": context, "url": url, "steps": []}

    failing = failing_steps(steps)
    if not failing:
        entry["error"] = "no step reported a nonzero return code"
        return entry

    for name, code, text in failing:
        errors = extract_step_errors(text, max_diagnostics)
        step = {
            "name": name,
            "exit_code": code,
            "kind": errors["kind"],
            "killed_after": errors["killed_after"],
        }

        container = extract_container(text)
        step["container"] = container

        if errors["kind"] == "check":
            signatures = errors["signatures"]
            step["tests_failed"] = len(errors["failures"])
            step["signatures_total"] = len(signatures)
            step["signatures"] = [
                {"message": message, "count": count}
                for message, count in signatures.most_common(max_signatures)
            ]
            hint_texts.update(signatures)

            for failure in errors["failures"]:
                test_entry = rollup.setdefault(
                    failure["test"],
                    {"status": failure["status"], "failures": [], "reasons": set()},
                )
                if failure["reason"]:
                    test_entry["reasons"].add(failure["reason"])
                test_entry["failures"].append(
                    Failure(
                        job=entry["context"],
                        command=failure["command"],
                        mode=failure["mode"] or "",
                        container=container,
                    )
                )
        else:
            step["diagnostics"] = errors["items"]
            step["diagnostics_total"] = errors["total"]
            hint_texts.update(errors["all_items"])
            if errors["kind"] == "hang":
                hint_texts.add(HANG)

        entry["steps"].append(step)

    return entry


def collect_errors(
    failed_jobs: Sequence[dict],
    max_jobs: int,
    max_diagnostics: int,
    max_signatures: int,
) -> dict:
    """
    Read the logs of each failed job and gather what they contain.

    Collection is kept separate from reporting so that the text digest and the
    JSON payload describe exactly the same findings.
    """
    jobs: List[dict] = []
    rollup: dict = {}
    hint_texts: set = set()

    for status in failed_jobs[:max_jobs]:
        context = status.get("context") or "?"
        url = status.get("target_url") or ""
        try:
            steps = fetch_job_steps(url)
        except Exception as e:  # noqa: BLE001 - a log we cannot read is not fatal
            jobs.append(
                {
                    "context": context,
                    "url": url,
                    "steps": [],
                    "error": f"could not read logs: {e}",
                }
            )
            continue

        jobs.append(
            collect_job_errors(
                context, url, steps, rollup, hint_texts, max_diagnostics, max_signatures
            )
        )

    return {
        "jobs": jobs,
        "rollup": rollup,
        "hints": hints_for(sorted(hint_texts)),
        "jobs_not_read": max(0, len(failed_jobs) - max_jobs),
    }


def collect_one_job(
    job_url: str,
    steps: Sequence[tuple],
    env: dict,
    max_diagnostics: int,
    max_signatures: int,
) -> dict:
    """
    Gather the findings of a single job, in the shape print_errors expects.

    Reporting on one job is the case where CIVET is reached by a job URL, which
    carries no pull request or commit to look the event up by.
    """
    rollup: dict = {}
    hint_texts: set = set()
    entry = collect_job_errors(
        env.get("CIVET_RECIPE_NAME") or "?",
        job_url,
        steps,
        rollup,
        hint_texts,
        max_diagnostics,
        max_signatures,
    )
    return {
        "jobs": [entry],
        "rollup": rollup,
        "hints": hints_for(sorted(hint_texts)),
        "jobs_not_read": 0,
    }



def print_errors(collected: dict, max_failures: int, label_jobs: bool = True) -> None:
    """
    Print the collected build diagnostics, error clusters and test rollup.

    Pass label_jobs=False when the caller has already named the job, as the
    single-job report does when it prints the event the job belongs to.
    """
    for job in collected["jobs"]:
        if label_jobs:
            print(f"\n{job['context']}  {job['url']}")
        if error := job.get("error"):
            print(f"  {error}")

        for step in job["steps"]:
            notes = [f"exit {step['exit_code']}"]
            if step["killed_after"] is not None:
                notes.append(f"killed by CIVET after {step['killed_after']} s")
            if step["kind"] == "check":
                notes.append(f"{step['tests_failed']} failing test(s) or example(s)")
                notes.append(f"{step['signatures_total']} error signature(s)")
                print(f"  {step['name']} ({', '.join(notes)})")
                for signature in step["signatures"]:
                    print(f"    x{signature['count']}  {signature['message'][:150]}")
                hidden = step["signatures_total"] - len(step["signatures"])
                if hidden:
                    print(f"    ... {hidden} more signature(s)")
            else:
                notes.append(step["kind"])
                print(f"  {step['name']} ({', '.join(notes)})")
                for item in step["diagnostics"]:
                    print(f"    {item[:200]}")
                hidden = step["diagnostics_total"] - len(step["diagnostics"])
                if hidden:
                    print(f"    ... {hidden} more unique")

    rollup = collected["rollup"]
    if rollup:
        print(f"\nunique failures: {len(rollup)}")
        for entry in rollup_entries(rollup, max_failures):
            print(f"  {entry['test']}  [{entry['status']}]")
            for reason in entry["reasons"][:3]:
                print(f"    {reason[:200]}")
            scope = (
                f"    failed in {len(entry['jobs'])} job(s), {entry['modes']} mode(s)"
            )
            if len(entry["jobs"]) <= 3:
                scope += f": {', '.join(entry['jobs'])}"
            print(scope)
            print(f"    {entry['reproduce']}")
            if entry["container"]:
                print(f"    in {entry['container']}")
            for hint in entry["hints"]:
                print(f"    -> {hint.split(' -> ', 1)[-1]}")
        if len(rollup) > max_failures:
            print(f"  ... {len(rollup) - max_failures} more; raise --max-failures")

    if collected["hints"]:
        print("\nknown remediations (job level):")
        for hint in collected["hints"]:
            print(f"  {hint}")

    if collected["jobs_not_read"]:
        print(
            f"\n... {collected['jobs_not_read']} more failed jobs; "
            "raise --max-log-jobs"
        )


def print_step_list(job_url: str, steps: Sequence[tuple]) -> None:
    """
    List a job's steps and their sizes.

    Sizes are what say whether a step can be printed at all, so they come
    before any request for its contents.
    """
    print(f"steps in {job_url} (pass --step NAME to read one):")
    for name, text in steps:
        print(f"  {name}  ({len(text) / 1024:.0f} KB)")


def grepped_lines(body: Sequence[str], pattern: str, context: int) -> List[int]:
    """Get the indices of the lines matching a pattern, plus their context."""
    regex = re.compile(pattern)
    keep: set = set()
    for i, line in enumerate(body):
        if regex.search(line):
            keep.update(range(max(0, i - context), min(len(body), i + context + 1)))
    return sorted(keep)


def print_step_log(
    steps: Sequence[tuple],
    step: str,
    lines: int,
    pattern: Optional[str],
    context: int,
) -> None:
    """
    Print the part of one step's log that was asked for, bounded to 'lines'.

    A tail is the wrong tool for a failure in a check step: `make check` keeps
    running the remaining unit tests and examples after one fails, so the
    detail sits in the middle of the log and no affordable tail reaches it.
    That is what the pattern is for. Either way at most 'lines' lines are
    printed, because a whole step can be several megabytes.
    """
    matches = [(n, t) for n, t in steps if step in n]
    if not matches:
        available = ", ".join(n for n, _ in steps)
        raise GitHubError(f"No step matching '{step}'. Available: {available}")

    for name, text in matches:
        body = text.splitlines()

        if not pattern:
            tail = body[-lines:]
            print(f"--- {name} (last {len(tail)} of {len(body)} lines) ---")
            print("\n".join(tail))
            continue

        keep = grepped_lines(body, pattern, context)
        shown = keep[:lines]
        print(
            f"--- {name} ({len(keep)} of {len(body)} lines match "
            f"'{pattern}', showing {len(shown)}) ---"
        )
        previous = None
        for i in shown:
            # A gap in the kept indices is a jump in the log, and hiding it
            # would make unrelated lines read as consecutive
            if previous is not None and i != previous + 1:
                print("...")
            print(body[i])
            previous = i
        if len(keep) > len(shown):
            print(f"... {len(keep) - len(shown)} more matching line(s); raise --lines")


def load_recipe_graph(recipes_dir: str) -> Optional[dict]:
    """
    Build the recipe dependency graph from a civet_recipes checkout.

    Entirely optional: CIVET's job statuses say nothing about dependencies, so
    a pending job that can never run looks identical to one that is merely
    queued. The recipes are the only place that relation is recorded, and not
    everyone has access to them.
    """
    files = glob.glob(os.path.join(recipes_dir, "**", "*.cfg"), recursive=True)
    if not files:
        print(
            f"WARNING: no recipes found under {recipes_dir}; "
            "dependency analysis skipped",
            file=sys.stderr,
        )
        return None

    files_of: dict = collections.defaultdict(set)
    deps: dict = {}
    for path in files:
        parser = configparser.ConfigParser(strict=False, interpolation=None)
        try:
            parser.read(path)
        except configparser.Error:
            continue
        if not parser.has_section("Main"):
            continue

        main = parser["Main"]
        name = (main.get("display_name") or main.get("name") or "").strip()
        if not name:
            continue

        # Dependency paths in a recipe are relative to the checkout root
        key = os.path.relpath(path, recipes_dir)
        files_of[name].add(key)
        deps[key] = (
            [
                value.strip()
                for option, value in parser["PullRequest Dependencies"].items()
                if option.startswith("filename")
            ]
            if parser.has_section("PullRequest Dependencies")
            else []
        )

    return {"files_of": dict(files_of), "deps": deps}


def downstream_of(graph: dict, contexts: Sequence[str]) -> set:
    """
    Get the display names that depend, directly or transitively, on the given ones.

    A display name can appear in more than one branch of the recipes, so a name
    counts as downstream if any of its recipe files is. That errs toward
    reporting a job as blocked rather than missing one.
    """
    dependents: dict = collections.defaultdict(set)
    for path, requires in graph["deps"].items():
        for required in requires:
            dependents[required].add(path)

    seen = set()
    frontier = {f for name in contexts for f in graph["files_of"].get(name, ())}
    while frontier:
        following = set()
        for path in frontier:
            for dependent in dependents.get(path, ()):
                if dependent not in seen:
                    seen.add(dependent)
                    following.add(dependent)
        frontier = following

    return {name for name, paths in graph["files_of"].items() if paths & seen}


def blocked_pending(graph: Optional[dict], statuses: Sequence[dict]) -> List[str]:
    """Get the pending jobs that are downstream of a failure, if recipes allow."""
    if not graph:
        return []
    failed = [s.get("context") for s in real_failures(statuses)]
    blocked = downstream_of(graph, [c for c in failed if c])
    return sorted(
        {
            s.get("context")
            for s in pending_jobs(statuses)
            if s.get("context") in blocked
        }
    )


def pending_jobs(statuses: Sequence[dict]) -> List[dict]:
    """
    Select the jobs that have not finished.

    Completeness is judged from the statuses that have been posted, so a job
    CIVET has not reported on at all is invisible here.
    """
    return [s for s in statuses if s.get("state") == PENDING_STATE]


def real_failures(statuses: Sequence[dict]) -> List[dict]:
    """
    Select the jobs that actually failed.

    Jobs CIVET skipped because a dependency failed also report a failure
    state, but there is nothing in them to fix.
    """
    return [
        s
        for s in statuses
        if s.get("state") in ("failure", "error")
        and s.get("description") != BLOCKED_DESCRIPTION
    ]


def describe_target(info: dict) -> str:
    """Name what is being reported on, which may be a commit rather than a PR."""
    sha = info["headRefOid"][:7]
    if info.get("number") is None:
        return f"commit {sha}"
    return f"PR #{info['number']} {info['headRefName']} @ {sha}"


def print_digest(info: dict, state: dict, graph: Optional[dict] = None) -> None:
    """Print the job-level report."""
    statuses = state["statuses"]
    failed_jobs = real_failures(statuses)
    blocked = [s for s in statuses if s.get("description") == BLOCKED_DESCRIPTION]

    pending = pending_jobs(statuses)
    passed = [s for s in statuses if s.get("state") == "success"]

    counts = [f"{len(failed_jobs)} failed"]
    if blocked:
        counts.append(f"{len(blocked)} blocked")
    if pending:
        counts.append(f"{len(pending)} pending")
    counts.append(f"{len(passed)} passed")
    print(f"{describe_target(info)}  jobs: {', '.join(counts)} / {len(statuses)}")

    if pending:
        # Every conclusion below is provisional while jobs are still running,
        # and GitHub's own rollup state does not reveal this
        print(
            f"\nEVENT INCOMPLETE: {len(pending)} job(s) still running. "
            "Failures seen so far are not the whole picture."
        )

    doomed = blocked_pending(graph, statuses)
    if doomed:
        print(
            f"\n{len(doomed)} of the {len(pending)} pending job(s) depend on a "
            "failed job and will not run:"
        )
        for name in doomed[:MAX_JOBS]:
            print(f"  {name}")
        if len(doomed) > MAX_JOBS:
            print(f"  ... and {len(doomed) - MAX_JOBS} more")

    if failed_jobs:
        print("\nfailed jobs:")
        for status in failed_jobs[:MAX_JOBS]:
            print(f"  {status.get('context')}  {status.get('target_url')}")
        if len(failed_jobs) > MAX_JOBS:
            print(f"  ... and {len(failed_jobs) - MAX_JOBS} more")

    if blocked:
        # Cascade effects of the failures above, with nothing of their own to fix
        print(
            "\nblocked by failed dependencies: "
            + ", ".join(s.get("context") or "?" for s in blocked)
        )


def json_payload(
    info: dict,
    state: dict,
    collected: dict,
    max_failures: int,
    graph: Optional[dict] = None,
) -> dict:
    """
    Build the machine-readable form of the same report.

    Capped the same way as the text output, with exact totals alongside, so
    that a consumer sees a bounded payload and still knows what was left out.
    """
    return {
        "pr": info["number"],
        "branch": info["headRefName"],
        "head_sha": info["headRefOid"],
        # GitHub's rollup, which reads "failure" while jobs are still pending
        "github_state": state["state"],
        "complete": bool(state["statuses"]) and not pending_jobs(state["statuses"]),
        "counts": {
            "jobs": len(state["statuses"]),
            "failed": len(real_failures(state["statuses"])),
            "blocked": sum(
                1
                for s in state["statuses"]
                if s.get("description") == BLOCKED_DESCRIPTION
            ),
            "pending": len(pending_jobs(state["statuses"])),
            "pending_blocked": len(blocked_pending(graph, state["statuses"])),
            "passed": sum(1 for s in state["statuses"] if s.get("state") == "success"),
        },
        # Passing jobs are omitted; they are counted above and carry nothing
        # a consumer would act on
        "jobs": [
            {
                "context": s.get("context"),
                "state": s.get("state"),
                "description": s.get("description"),
                "url": s.get("target_url"),
                "blocked": s.get("description") == BLOCKED_DESCRIPTION,
            }
            for s in state["statuses"]
            if s.get("state") != "success"
        ],
        "failed_jobs": collected["jobs"],
        "tests_total": len(collected["rollup"]),
        "tests": rollup_entries(collected["rollup"], max_failures),
        "hints": collected["hints"],
        "jobs_not_read": collected["jobs_not_read"],
        # Empty unless a recipes checkout was given with --recipes
        "pending_blocked": blocked_pending(graph, state["statuses"]),
    }


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    target = parser.add_mutually_exclusive_group()
    target.add_argument("--pr", type=int, help="PR number (default: branch's PR)")
    target.add_argument(
        "--sha",
        type=str,
        help="Report on a commit instead of a pull request, which is what makes "
        "pushes to devel reportable",
    )
    target.add_argument(
        "--job",
        type=str,
        metavar="JOB_URL",
        help="Report on one CIVET job named by its URL. Reads that job's logs "
        "and reports which event and commit it belongs to",
    )
    parser.add_argument(
        "--repo",
        type=str,
        help="owner/name of the repository (default: the upstream remote's)",
    )
    parser.add_argument(
        "--upstream-remote",
        type=str,
        help="Remote naming the upstream repo. Defaults to CIVET_UPSTREAM_REMOTE, "
        f"else the first of {', '.join(UPSTREAM_REMOTE_CANDIDATES)} that exists",
    )
    parser.add_argument(
        "--max-failures",
        type=int,
        default=DEFAULT_MAX_FAILURES,
        help=f"Max failures to describe (default: {DEFAULT_MAX_FAILURES})",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        dest="as_json",
        help="Emit the report as JSON; implies --errors",
    )
    parser.add_argument(
        "--recipes",
        type=str,
        metavar="DIR",
        help="Optional civet_recipes checkout, used to report which pending "
        "jobs depend on a failed job and so will never run",
    )
    parser.add_argument(
        "--errors",
        action="store_true",
        help="Pull the logs of each failed job and extract build/test errors",
    )
    parser.add_argument(
        "--max-log-jobs",
        type=int,
        default=DEFAULT_MAX_LOG_JOBS,
        help=f"With --errors, max jobs to read (default: {DEFAULT_MAX_LOG_JOBS})",
    )
    parser.add_argument(
        "--max-diagnostics",
        type=int,
        default=DEFAULT_MAX_DIAGNOSTICS,
        help=f"With --errors, unique errors per step (default: {DEFAULT_MAX_DIAGNOSTICS})",
    )
    parser.add_argument(
        "--max-signatures",
        type=int,
        default=DEFAULT_MAX_SIGNATURES,
        help=f"With --errors, error signatures per step (default: {DEFAULT_MAX_SIGNATURES})",
    )
    parser.add_argument(
        "--list-steps",
        action="store_true",
        help="With --job, list the job's steps and their sizes",
    )
    parser.add_argument(
        "--step",
        type=str,
        help="With --job, print the log of the step whose name contains this",
    )
    parser.add_argument(
        "--grep",
        type=str,
        metavar="PATTERN",
        help="With --step, print the lines matching this regular expression "
        "instead of a tail, which is the only way to reach a failure in the "
        "middle of a large log",
    )
    parser.add_argument(
        "--context",
        type=int,
        default=0,
        help="With --grep, lines of context to print around each match, 0 by default",
    )
    parser.add_argument(
        "--lines",
        type=int,
        default=DEFAULT_MAX_LOG_LINES,
        help="With --step, the most lines to print, either as a tail or as "
        f"matches (default: {DEFAULT_MAX_LOG_LINES})",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    """Perform the main action; run from __main__."""
    args = parse_args(argv)

    if args.job:
        # Fetched once and reused, because a job's whole log is a single
        # tarball download regardless of which part of it is wanted
        try:
            steps = fetch_job_steps(args.job)
        except Exception as e:  # noqa: BLE001 - reported, not raised as a trace
            raise GitHubError(f"could not read the logs of {args.job}: {e}") from e

        if args.list_steps:
            print_step_list(args.job, steps)
        elif args.step:
            print_step_log(steps, args.step, args.lines, args.grep, args.context)
        else:
            env = extract_job_info(steps)
            print_job_header(args.job, env)
            print_errors(
                collect_one_job(
                    args.job, steps, env, args.max_diagnostics, args.max_signatures
                ),
                args.max_failures,
                label_jobs=False,
            )
        return 0

    slug = resolve_repo(args.repo, args.upstream_remote)
    if not slug:
        tried = args.upstream_remote or os.environ.get("CIVET_UPSTREAM_REMOTE")
        tried = [tried] if tried else list(UPSTREAM_REMOTE_CANDIDATES)
        raise GitHubError(
            f"No repository found from remote(s) {', '.join(tried)}; "
            "pass --repo owner/name or set --upstream-remote"
        )

    info = resolve_target(slug, args.pr, args.sha)
    sha = info["headRefOid"]

    state = fetch_state(slug, sha)

    # Optional: without it, a doomed pending job is indistinguishable from a
    # queued one, which is a limit of the status data rather than a bug
    graph = load_recipe_graph(args.recipes) if args.recipes else None

    # JSON is for programmatic use, where the per-test detail is the point
    collect = args.errors or args.as_json
    collected = (
        collect_errors(
            real_failures(state["statuses"]),
            args.max_log_jobs,
            args.max_diagnostics,
            args.max_signatures,
        )
        if collect
        else None
    )

    if args.as_json:
        print(
            json.dumps(
                json_payload(info, state, collected, args.max_failures, graph),
                indent=1,
                sort_keys=True,
            )
        )
        return 0

    print_digest(info, state, graph)
    if collected:
        print_errors(collected, args.max_failures)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
