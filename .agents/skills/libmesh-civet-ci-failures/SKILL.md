---
name: libmesh-civet-ci-failures
description: >-
  Use when investigating why CIVET failed or timed out for a libMesh pull request, commit, or
  single CIVET job URL: which jobs failed, which step and process count, which unit test or
  example, and how to reproduce it locally. Reads GitHub commit statuses and CIVET job logs
  through the bundled scripts/civet_ci_failures.py, which understands `make check` output: CppUnit
  failures, failed examples, and hangs that CIVET kills at its step time limit.
---

# libMesh CIVET failures

## The tool

`scripts/civet_ci_failures.py`, in this skill's directory, needs an authenticated `gh` CLI and
network access to `civet.inl.gov`. It combines the commit statuses CIVET posts, one per job, with
the job step logs CIVET serves as a public tarball, which is where build errors and test failures
live. Run it from a libMesh checkout, whose remotes it uses to find the repository:

```bash
B=.agents/skills/libmesh-civet-ci-failures/scripts/civet_ci_failures.py
$B                                   # the current branch's PR, job level only
$B --pr 4574 --errors                # a PR, reading each failed job's logs
$B --pr 4574 --json                  # the same findings as JSON; implies --errors
$B --sha <commit> --errors           # a commit, e.g. a push to devel
$B --job https://civet.inl.gov/job/4223219/   # one job: recipe, PR, commits, findings
```

The repository comes from `--repo owner/name`, else the remote named by `--upstream-remote` or
`$CIVET_UPSTREAM_REMOTE`, else the first of `up`, `upstream`, `origin` that exists. If none of
those points at libMesh/libmesh, pass `--repo libMesh/libmesh`.

`--job` reports one job in isolation, so it cannot tell whether other jobs failed the same way.
Widen to the event with `--pr` or `--sha` on the identifiers its header reports. A failure
reported by several jobs is merged into one entry with a count of jobs and of modes, where a mode
is a distinct invocation, such as a different process count.

The tool is originally adapted from MOOSE's `python/civet_ci_failures`; the job-level logic is shared, and
the log parsing is libMesh's own. Its unit tests are in `scripts/test_civet_ci_failures.py`:

```bash
python3 -m unittest discover -s .agents/skills/libmesh-civet-ci-failures/scripts -p 'test_civet_ci_failures.py'
```

## Reading a raw log

Only for what the report does not cover. Start with `--list-steps` for sizes, then read one step:

```bash
$B --job <url> --list-steps
$B --job <url> --step 08_Check_7 --lines 200          # the tail of one step
$B --job <url> --step 08_Check_7 --grep '^[0-9]+\) test: ' --context 3
```

`--grep` prints matches from the top of the log, capped by `--lines`, so a pattern that also
matches build output early in the step can use up the cap before reaching the tests. Use a pattern
specific to what you want, or `--lines` alone for the tail. To search one job's logs repeatedly,
download the tarball once:

```bash
mkdir -p /tmp/civet/<id> && curl -sSL https://civet.inl.gov/job_results/<id>/ | tar xz -C /tmp/civet/<id>
```

Never read a whole step into context; a check step is megabytes.

## Reading the job states correctly

- **A "failed" job is often a cascade.** CIVET marks a job skipped because a dependency failed
  with `Won't run due to failed dependencies`, and GitHub reports it as a failure. The report
  counts these separately as blocked.
- **GitHub's combined state is not a completion signal.** It turns `failure` as soon as one job
  fails, with the others still pending. The report says `EVENT INCOMPLETE` while any are pending.
- **Absence of pending statuses is not completion either.** A job that has not started has no
  status at all, so require the counts to hold still across two polls before believing an event
  finished.
- **Pushing to a PR cancels the in-flight event**, so batch fixes into one push unless an early
  failure has already doomed the event.

## What a libMesh check step runs

A `Check_N` step runs the full `make check` with `LIBMESH_RUN` set to `mpiexec -n N`: the bundled
contrib test suites (netCDF and others), then `tests/unit_tests-<method>`, then every example's
`run.sh`. The sweep recipes, `Distributed make check sweep (odd)` and `(even)`, run one such step
per process count in ascending order and stop at the first that fails, so the later steps have
empty logs. A failure appearing only in a sweep, at process counts the rest of the matrix never
uses, points to behavior that depends on the partition.

The report's reproduce commands are the invocations the step logged, with `cd` paths relative to
the build directory: `tests` for the unit tests, `examples/<category>/<name>` for an example.
They carry the PETSc options CI passes, which can matter.

## A hang looks like a crash

A hung step runs until CIVET's step limit, then CIVET kills it. The log then ends with CIVET's
cancel notice, `Caught signal number 15 Terminate` from PETSc's signal handler on each rank,
`MPI_Abort`, and `make` reporting `Terminated`, with exit code 143. These lines come from the
kill, so the report leaves them out of the error signatures and labels the step `killed by CIVET
after N s`.

The unit test driver prints test names only on failure. As each test starts CppUnit prints a dot,
so a hung run ends in a line of dots with no name, and the report names the hung test by position:
`unit_tests-opt test #9` means the 9th test the run started. To turn the position into a name,
run the reported command, which adds `--verbose` to print each test's name as it starts. The
position in the source order is a strong hint as well: CppUnit orders suites by the address of
their static registration objects, which in practice follows `unit_tests_SOURCES` in
`tests/Makefile.am`, and orders tests within a suite by `CPPUNIT_TEST`, minus those that
configure-dependent `#if` guards exclude. When the same position appears at several process
counts, that agreement is itself good evidence.

## Why a test hangs only at high process counts

A failed CppUnit assertion throws an exception. If the assertion fails on only some ranks, only
those ranks throw, and their stacks unwind out of the test method. The other ranks keep executing
the test and block in its next collective operation: a vector `close()`,
`System::update()`, `create_dof_constraints()`, a `comm().sum()`, anything marked
`parallel_object_only()`. Nothing prints, because CppUnit writes failures only when the run ends,
and only rank 0's output reaches the log.

The usual reason an assertion holds on every rank at low process counts and fails on one at high
counts is a small test mesh. As the process count grows, some rank owns no elements, or only
interior elements, or no elements on a given boundary, and a check like `CPPUNIT_ASSERT(!values.empty())`
on rank-local data fails there. Look in the identified test for assertions on rank-local
quantities that precede a collective call. Fix them by reducing over the communicator before
asserting, or by asserting a property that holds on every rank by construction.

The same rank-0-only output means a unit test can fail on another rank while rank 0 prints `OK`;
the run then exits nonzero with no failure in the log. Rerun locally with `--keep-cout` to see
every rank's output.

## Reproducing

Run the identified suite at the failing process count, selecting it with `--re` (a bare suite
name is not a filter):

```bash
cd <build>/tests && mpiexec -n <N> ./unit_tests-opt --re 'DofMapTest'
```

Reproduce a hang with a dbg build when one is available. CI's sweeps build opt, where ranks that
diverge simply wait on each other. In dbg the `parallel_object_only()` checks compare ranks, so the
same divergence aborts at once with ``Assertion `(this->comm()).verify(...)' failed``, naming the
rank and the file and line where the mismatch was detected. On the throwing rank that location is
usually a destructor run during stack unwinding, such as `PetscVector::clear()`, whose
`parallel_object_only()` check fails because the other ranks are still in the test body. The
defect is the failed assertion which only threw on some ranks earlier in the same test.

The sweep recipes build with a distributed mesh (`--enable-distmesh`, which the recipes still
spell with its deprecated alias `--enable-parmesh`), so a local build using a replicated mesh
can partition differently and may not reproduce a partition-dependent failure. Check the job's
`Configure` step for the configure line before concluding that a failure does not reproduce.

## Limits

- A failed contrib test is not parsed; the step falls back to the build-failure reading, which
  reports make's complaint. Grep the step for `FAIL: ` to find it.
- The unit tests cover the parsing against synthetic logs shaped like real ones; a change to
  CppUnit's, `run_common.sh`'s or CIVET's output format can still slip past them.

## When something does not fit

Ask the user rather than guessing. A recipe whose steps are not `make check`, or a log the report
misreads, means the tool or this skill is out of date: say what differs, and update both, with a
unit test for the new log shape, once the right reading is established.
