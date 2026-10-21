#!/usr/bin/env python3
"""
Unit tests for civet_ci_failures.py. Run with

  python3 -m unittest discover -s .agents/skills/libmesh-civet-ci-failures/scripts \
    -p 'test_civet_ci_failures.py'

Adapted from the tests of MOOSE's python/civet_ci_failures.
"""

import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile
import unittest
from contextlib import redirect_stdout
from unittest import mock

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

import civet_ci_failures as ccf


class TestLatestStatuses(unittest.TestCase):
    def testNewestPerContextWins(self):
        # GitHub returns every status ever posted, newest first; CIVET posts
        # one when a job starts and another when it finishes.
        statuses = [
            {"context": "build", "state": "success"},
            {"context": "test", "state": "pending"},
            {"context": "build", "state": "pending"},
        ]
        self.assertEqual(
            ccf.latest_statuses(statuses),
            [
                {"context": "build", "state": "success"},
                {"context": "test", "state": "pending"},
            ],
        )

    def testMissingContextKeptOnce(self):
        statuses = [{"state": "success"}, {"state": "pending"}]
        self.assertEqual(ccf.latest_statuses(statuses), [{"state": "success"}])


class TestPendingAndRealFailures(unittest.TestCase):
    def testPendingJobs(self):
        statuses = [
            {"context": "a", "state": "pending"},
            {"context": "b", "state": "success"},
        ]
        self.assertEqual(ccf.pending_jobs(statuses), [statuses[0]])

    def testRealFailuresExcludesCascades(self):
        statuses = [
            {"context": "a", "state": "failure", "description": "boom"},
            {
                "context": "b",
                "state": "failure",
                "description": ccf.BLOCKED_DESCRIPTION,
            },
            {"context": "c", "state": "error", "description": "boom"},
            {"context": "d", "state": "success", "description": None},
        ]
        self.assertEqual(ccf.real_failures(statuses), [statuses[0], statuses[2]])


class TestDescribeTarget(unittest.TestCase):
    def testPullRequest(self):
        info = {
            "number": 123,
            "headRefName": "my-branch",
            "headRefOid": "abcdef0123456",
        }
        self.assertEqual(ccf.describe_target(info), "PR #123 my-branch @ abcdef0")

    def testBareCommit(self):
        info = {"number": None, "headRefName": None, "headRefOid": "abcdef0123456"}
        self.assertEqual(ccf.describe_target(info), "commit abcdef0")


class TestNormalizeError(unittest.TestCase):
    def testTmpDirAddressAndDigitsAreNormalized(self):
        message = "opening /tmp/abc123/file at 0xdeadBEEF, attempt 3"
        self.assertEqual(
            ccf.normalize_error(message),
            "opening /tmp/<dir>/file at <addr>, attempt N",
        )



UNIT_INVOCATION = "mpiexec -n 7 ./unit_tests-opt --option-with-dashes -pc_type asm"

# A check step whose unit tests hung: the progress marks are interleaved with
# test output, and rank 0's last marks are flushed only as the kill arrives
HUNG_UNIT_TEXT = (
    "Executing step in docker://img:tag\n"
    "make[3]: Entering directory '/tmp/build/libmesh/build/tests'\n"
    "make[3]: warning: overriding recipe for target 'unit_tests-kokkos-dbg'\n"
    f"{UNIT_INVOCATION}\n"
    "Will run the following tests:\n"
    "All Tests\n"
    "....\n"
    "Warning!  Parallel xda/xdr is not yet implemented.\n"
    "..Writing a serialized file instead.\n"
    "\n"
    "*****************************************************\n"
    "CIVET: Cancelling job due to step taking longer than the max 21600 seconds\n"
    "\n"
    "...[0]PETSC ERROR: Caught signal number 15 Terminate: Some process has told this process to end\n"
    "application called MPI_Abort(MPI_COMM_WORLD, 59) - process 0\n"
    "make: *** [Makefile:36055: check-recursive] Terminated\n"
    "completed with return code 143\n"
)

# A check step whose unit tests finished with one failed assertion and one
# uncaught exception
FAILED_UNIT_TEXT = (
    f"{UNIT_INVOCATION}\n"
    "Will run the following tests:\n"
    "All Tests\n"
    "..F..E.\n"
    "\n"
    "!!!FAILURES!!!\n"
    "Test Results:\n"
    "Run:  5   Failures: 1   Errors: 1\n"
    "\n"
    "1) test: DofMapTest::testFoo (F) line: 437 ../tests/base/dof_map_test.C\n"
    "assertion failed\n"
    "- Expression: !values.empty()\n"
    "\n"
    "2) test: MeshTest::testBar (E)\n"
    "uncaught exception of type std::exception (or derived).\n"
    "- libMesh error\n"
    "\n"
    "make[3]: *** [Makefile:13855: check-TESTS] Error 1\n"
    "completed with return code 2\n"
)

EXAMPLE_COMMAND = "mpiexec -n 5 ./example-opt -d 3 mesh.xda"

# A check step in which one example ran to completion and the next died on a
# failed assertion
FAILED_EXAMPLE_TEXT = (
    f"{UNIT_INVOCATION}\n"
    ".....\n"
    "OK (5 tests)\n"
    "make[2]: Entering directory '/tmp/build/libmesh/build/examples/introduction/introduction_ex1'\n"
    "***************************************************************\n"
    "* Running Example introduction_ex1:\n"
    f"*  {EXAMPLE_COMMAND}\n"
    "***************************************************************\n"
    "* Done Running Example introduction_ex1:\n"
    f"*  {EXAMPLE_COMMAND}\n"
    "***************************************************************\n"
    "PASS: run.sh\n"
    "make[2]: Entering directory '/tmp/build/libmesh/build/examples/introduction/introduction_ex2'\n"
    "***************************************************************\n"
    "* Running Example introduction_ex2:\n"
    "*  mpiexec -n 5 ./example-opt\n"
    "***************************************************************\n"
    "Assertion `i < size()' failed.\n"
    "Assertion `i < size()' failed.\n"
    "FAIL: run.sh\n"
    "completed with return code 2\n"
)


class TestParseRuns(unittest.TestCase):
    def testUnitRunAndExamplesWithTheirCommandsAndDirectories(self):
        runs = ccf.parse_runs(FAILED_EXAMPLE_TEXT.splitlines())
        self.assertEqual(
            [(r["kind"], r["name"], r["finished"]) for r in runs],
            [
                ("unit", "unit_tests-opt", True),
                ("example", "introduction_ex1", True),
                ("example", "introduction_ex2", False),
            ],
        )
        self.assertEqual(runs[0]["command"], UNIT_INVOCATION)
        self.assertEqual(runs[0]["directory"], "tests")
        self.assertEqual(runs[1]["command"], EXAMPLE_COMMAND)
        self.assertEqual(runs[1]["directory"], "examples/introduction/introduction_ex1")

    def testSerialInvocationWithoutLauncher(self):
        runs = ccf.parse_runs(["./unit_tests-dbg --re Foo", "OK (1 test)"])
        self.assertEqual(runs[0]["command"], "./unit_tests-dbg --re Foo")
        self.assertTrue(runs[0]["finished"])

    def testBuildingTheBinaryIsNotARun(self):
        lines = [
            "make  unit_tests-opt",
            "libtool: link: mpicxx -O2 -o .libs/unit_tests-opt driver.o",
        ]
        self.assertEqual(ccf.parse_runs(lines), [])


class TestProgressMarks(unittest.TestCase):
    def testCountsMarksAtLineStartsAcrossInterleavedOutput(self):
        lines = ["....", "Warning!", "..Writing a file", "...[0]PETSC ERROR: x"]
        self.assertEqual(ccf.progress_marks(lines), (9, 0))

    def testFailureMarksCountedButWordsAreNot(self):
        lines = ["..F..E.", "Failures: 1", "ERROR: x", "..Error in foo"]
        self.assertEqual(ccf.progress_marks(lines), (7, 2))


class TestCppunitFailures(unittest.TestCase):
    def testFailureAndErrorDetailsWithoutBoilerplate(self):
        failures = ccf.cppunit_failures(FAILED_UNIT_TEXT.splitlines())
        self.assertEqual(
            failures,
            [
                {
                    "test": "DofMapTest::testFoo",
                    "status": ccf.FAILURE,
                    "detail": "- Expression: !values.empty()",
                    "location": "../tests/base/dof_map_test.C:437",
                },
                {
                    "test": "MeshTest::testBar",
                    "status": ccf.ERROR,
                    "detail": "uncaught exception of type std::exception (or derived)."
                    " / - libMesh error",
                    "location": None,
                },
            ],
        )


class TestErrorSignatures(unittest.TestCase):
    def testRepeatedSignatureIsCountedAndPetscBoilerplateIsNot(self):
        lines = [
            "Assertion `i < size()' failed.",
            "Assertion `i < size()' failed.",
            "[0]PETSC ERROR: Try option -start_in_debugger",
            "[0]PETSC ERROR: Caught signal number 11 SEGV: Segmentation Violation",
        ]
        signatures = ccf.error_signatures(lines)
        self.assertEqual(signatures["Assertion `i < size()' failed."], 2)
        self.assertEqual(len(signatures), 2)


class TestExtractCheckStepErrors(unittest.TestCase):
    def testHungUnitTestIsNamedByItsPosition(self):
        result = ccf.extract_step_errors(HUNG_UNIT_TEXT, limit=5)
        self.assertEqual(result["kind"], "check")
        self.assertEqual(result["killed_after"], 21600)
        (failure,) = result["failures"]
        self.assertEqual(failure["test"], "unit_tests-opt test #9")
        self.assertEqual(failure["status"], ccf.HANG)
        self.assertEqual(failure["command"], f"cd tests && {UNIT_INVOCATION} --verbose")
        # The SIGTERM the kill delivers is not the cause
        self.assertEqual(len(result["signatures"]), 0)

    def testFailedUnitTestsCarryTheirDetailAndSelector(self):
        result = ccf.extract_step_errors(FAILED_UNIT_TEXT, limit=5)
        self.assertIsNone(result["killed_after"])
        self.assertEqual(
            [(f["test"], f["status"]) for f in result["failures"]],
            [("DofMapTest::testFoo", ccf.FAILURE), ("MeshTest::testBar", ccf.ERROR)],
        )
        self.assertEqual(
            result["failures"][0]["command"],
            f"cd tests && {UNIT_INVOCATION} --re '^DofMapTest::testFoo$'",
        )
        self.assertIn("dof_map_test.C:437", result["failures"][0]["reason"])

    def testFailedExampleCarriesItsSignatureAndDirectory(self):
        result = ccf.extract_step_errors(FAILED_EXAMPLE_TEXT, limit=5)
        (failure,) = result["failures"]
        self.assertEqual(failure["test"], "example introduction_ex2")
        self.assertEqual(failure["status"], ccf.FAILURE)
        self.assertEqual(failure["reason"], "Assertion `i < size()' failed.")
        self.assertEqual(
            failure["command"],
            "cd examples/introduction/introduction_ex2 && mpiexec -n 5 ./example-opt",
        )

    def testHungExampleIsReportedAsAHang(self):
        text = FAILED_EXAMPLE_TEXT.replace(
            "FAIL: run.sh\n",
            "CIVET: Cancelling job due to step taking longer than the max 600 seconds\n",
        )
        (failure,) = ccf.extract_step_errors(text, limit=5)["failures"]
        self.assertEqual(failure["test"], "example introduction_ex2")
        self.assertEqual(failure["status"], ccf.HANG)

    def testCrashedUnitRunWithoutResult(self):
        text = (
            f"{UNIT_INVOCATION}\n"
            "...\n"
            "[0]PETSC ERROR: Caught signal number 11 SEGV: Segmentation Violation\n"
            "completed with return code 2\n"
        )
        (failure,) = ccf.extract_step_errors(text, limit=5)["failures"]
        self.assertEqual(failure["test"], "unit_tests-opt test #3")
        self.assertEqual(failure["status"], ccf.CRASH)

    def testHangOutsideAnyRunFallsBackToTheTailBeforeTheKill(self):
        text = (
            "compiling foo.C\n"
            "CIVET: Cancelling job due to step taking longer than the max 600 seconds\n"
            "make: *** [all] Terminated\n"
            "completed with return code 143\n"
        )
        result = ccf.extract_step_errors(text, limit=5)
        self.assertEqual(result["kind"], "hang")
        self.assertEqual(result["items"], ["compiling foo.C"])


class TestFailingSteps(unittest.TestCase):
    def testOnlyNonzeroReturnCodesSelected(self):
        steps = [
            ("01_Build", "...\ncompleted with return code 0\n"),
            ("02_Test", "...\ncompleted with return code 1\n"),
        ]
        self.assertEqual(ccf.failing_steps(steps), [("02_Test", 1, steps[1][1])])

    def testLastReturnCodeWinsWhenStepRetried(self):
        text = (
            "completed with return code 1\n...retry...\ncompleted with return code 85\n"
        )
        self.assertEqual(ccf.failing_steps([("step", text)]), [("step", 85, text)])

    def testNoReturnCodeMeansNotFailing(self):
        self.assertEqual(ccf.failing_steps([("step", "no marker here")]), [])


class TestExtractStepErrors(unittest.TestCase):
    def testBuildKindFromCompilerDiagnostic(self):
        text = "foo.cc:12:3: error: use of undeclared identifier 'x'\n"
        result = ccf.extract_step_errors(text, limit=5)
        self.assertEqual(result["kind"], "build")
        self.assertEqual(result["total"], 1)

    def testBuildDiagnosticPreemptsMakeFailureSymptom(self):
        text = (
            "foo.cc:12:3: error: use of undeclared identifier 'x'\n"
            "make[1]: *** [foo.o] Error 1\n"
        )
        result = ccf.extract_step_errors(text, limit=5)
        self.assertEqual(result["kind"], "build")
        self.assertEqual(
            result["items"], ["foo.cc:12:3: error: use of undeclared identifier 'x'"]
        )

    def testInfraKindWinsOverNonDiagnosticBuildComplaint(self):
        text = (
            "fatal: unable to access 'https://github.com/libMesh/TIMPI.git/': "
            "error: RPC failed\n"
            "curl: (28) connection timed out\n"
        )
        result = ccf.extract_step_errors(text, limit=5)
        self.assertEqual(result["kind"], "infra")

    def testMakeFailureFallsBackToBuildWhenNothingElseMatches(self):
        text = "make[2]: *** [target] Error 2\n"
        result = ccf.extract_step_errors(text, limit=5)
        self.assertEqual(result["kind"], "build")

    def testUnknownKindFallsBackToTailWithoutTeardown(self):
        text = "useful line one\nuseful line two\nRemoving container\nSubmitting step statistics\n"
        result = ccf.extract_step_errors(text, limit=5)
        self.assertEqual(result["kind"], "unknown")
        self.assertEqual(result["items"], ["useful line one", "useful line two"])

    def testLimitCapsItemsButNotTotal(self):
        lines = "\n".join(f"error: problem {i}" for i in range(10))
        result = ccf.extract_step_errors(lines, limit=3)
        self.assertEqual(len(result["items"]), 3)
        self.assertEqual(result["total"], 10)


class TestRollupEntries(unittest.TestCase):
    def testShortestCommandChosenWithMatchingContainer(self):
        rollup = {
            "a.test": {
                "status": "ERROR",
                "reasons": {"TIMEOUT"},
                "failures": [
                    ccf.Failure(
                        job="job1",
                        command="cd tests && mpiexec -n 4 ./unit_tests-opt --re a",
                        mode="-p 4",
                        container="docker://x",
                    ),
                    ccf.Failure(
                        job="job2",
                        command="./unit_tests-opt --re a",
                        mode="-p 8",
                        container=None,
                    ),
                ],
            },
        }
        entries = ccf.rollup_entries(rollup, limit=10)
        self.assertEqual(len(entries), 1)
        entry = entries[0]
        self.assertEqual(entry["test"], "a.test")
        self.assertEqual(entry["jobs"], ["job1", "job2"])
        self.assertEqual(entry["modes"], 2)
        self.assertEqual(entry["reproduce"], "./unit_tests-opt --re a")
        self.assertIsNone(entry["container"])

    def testLimitCapsEntries(self):
        failure = ccf.Failure(
            job="job1", command="./unit_tests-opt --re x", mode="", container=None
        )
        rollup = {
            f"t{i}.test": {"status": "ERROR", "reasons": set(), "failures": [failure]}
            for i in range(5)
        }
        self.assertEqual(len(ccf.rollup_entries(rollup, limit=2)), 2)


class TestGreppedLines(unittest.TestCase):
    def testUnionOfMatchesWithContext(self):
        body = [
            "one",
            "two",
            "MATCH here",
            "four",
            "five",
            "six",
            "MATCH again",
            "eight",
        ]
        self.assertEqual(ccf.grepped_lines(body, "MATCH", 1), [1, 2, 3, 5, 6, 7])

    def testNoMatchIsEmpty(self):
        self.assertEqual(ccf.grepped_lines(["a", "b"], "nope", 0), [])


def cfg_file(directory, name, display_name, requires=()):
    """Write a minimal CIVET recipe .cfg with an optional dependency list."""
    lines = ["[Main]", f"display_name = {display_name}", ""]
    if requires:
        lines += ["[PullRequest Dependencies]"]
        lines += [f"filename{i + 1} = {r}" for i, r in enumerate(requires)]
    with open(os.path.join(directory, name), "w") as handle:
        handle.write("\n".join(lines) + "\n")


class TestLoadRecipeGraphAndDownstream(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        cfg_file(self.tmp, "build.cfg", "Build")
        cfg_file(self.tmp, "test.cfg", "Test", requires=["build.cfg"])
        cfg_file(self.tmp, "other.cfg", "Other")

    def testGraphMapsNamesAndDeps(self):
        graph = ccf.load_recipe_graph(self.tmp)
        self.assertEqual(
            graph["files_of"],
            {
                "Build": {"build.cfg"},
                "Test": {"test.cfg"},
                "Other": {"other.cfg"},
            },
        )
        self.assertEqual(graph["deps"]["test.cfg"], ["build.cfg"])

    def testDownstreamIsTransitiveDependent(self):
        graph = ccf.load_recipe_graph(self.tmp)
        self.assertEqual(ccf.downstream_of(graph, ["Build"]), {"Test"})

    def testBlockedPendingNeedsFailureAndGraph(self):
        graph = ccf.load_recipe_graph(self.tmp)
        statuses = [
            {"context": "Build", "state": "failure", "description": "x"},
            {"context": "Test", "state": "pending"},
            {"context": "Other", "state": "pending"},
        ]
        self.assertEqual(ccf.blocked_pending(graph, statuses), ["Test"])
        self.assertEqual(ccf.blocked_pending(None, statuses), [])

    def testEmptyDirWarnsAndReturnsNone(self):
        empty = tempfile.mkdtemp()
        with mock.patch("sys.stderr", new_callable=io.StringIO):
            self.assertIsNone(ccf.load_recipe_graph(empty))


TEST_STEP_TEXT = "Executing step in docker://img:tag\n" + FAILED_UNIT_TEXT


class TestCollectJobErrors(unittest.TestCase):
    def testTestStepPopulatesStepsAndRollup(self):
        steps = [
            ("01_Build", "completed with return code 0\n"),
            ("02_Test", TEST_STEP_TEXT),
        ]
        rollup, hint_texts = {}, set()
        entry = ccf.collect_job_errors(
            "ctx",
            "http://job/1",
            steps,
            rollup,
            hint_texts,
            max_diagnostics=5,
            max_signatures=5,
        )
        self.assertEqual(len(entry["steps"]), 1)
        step = entry["steps"][0]
        self.assertEqual(step["kind"], "check")
        self.assertEqual(step["tests_failed"], 2)
        self.assertEqual(step["container"], "docker://img:tag")
        self.assertIn("DofMapTest::testFoo", rollup)
        failure = rollup["DofMapTest::testFoo"]["failures"][0]
        self.assertEqual(
            failure.command,
            f"cd tests && {UNIT_INVOCATION} --re '^DofMapTest::testFoo$'",
        )
        self.assertEqual(failure.mode, UNIT_INVOCATION)

    def testNoFailingStepReportsError(self):
        entry = ccf.collect_job_errors(
            "ctx", "url", [("step", "completed with return code 0\n")], {}, set(), 5, 5
        )
        self.assertIn("error", entry)


class TestCollectOneJob(unittest.TestCase):
    def testUsesRecipeNameFromEnv(self):
        steps = [("01_Build", "completed with return code 1\nerror: something broke\n")]
        result = ccf.collect_one_job(
            "http://job/2", steps, {"CIVET_RECIPE_NAME": "My Recipe"}, 5, 5
        )
        self.assertEqual(result["jobs"][0]["context"], "My Recipe")
        self.assertEqual(
            result["jobs"][0]["steps"][0]["diagnostics"], ["error: something broke"]
        )
        self.assertEqual(result["jobs_not_read"], 0)


class TestCollectErrors(unittest.TestCase):
    def testFetchesEachFailedJobUpToMaxLogJobs(self):
        failed_jobs = [
            {"context": f"job{i}", "target_url": f"http://job/{i}"} for i in range(3)
        ]
        with mock.patch.object(
            ccf, "fetch_job_steps", return_value=[("step", TEST_STEP_TEXT)]
        ):
            collected = ccf.collect_errors(
                failed_jobs, max_jobs=2, max_diagnostics=5, max_signatures=5
            )
        self.assertEqual(len(collected["jobs"]), 2)
        self.assertEqual(collected["jobs_not_read"], 1)
        self.assertIn("MeshTest::testBar", collected["rollup"])

    def testUnreadableLogIsReportedNotRaised(self):
        failed_jobs = [{"context": "job0", "target_url": "http://job/0"}]
        with mock.patch.object(ccf, "fetch_job_steps", side_effect=RuntimeError("403")):
            collected = ccf.collect_errors(
                failed_jobs, max_jobs=5, max_diagnostics=5, max_signatures=5
            )
        self.assertIn("could not read logs", collected["jobs"][0]["error"])


SAMPLE_INFO = {"number": 42, "headRefName": "br", "headRefOid": "abc123def456"}
SAMPLE_STATE = {
    "state": "failure",
    "statuses": [
        {"context": "A", "state": "failure", "description": "boom", "target_url": "u1"},
        {"context": "B", "state": "success", "description": None, "target_url": "u2"},
        {"context": "C", "state": "pending", "description": None, "target_url": "u3"},
    ],
}


def captured(func, *args, **kwargs):
    buf = io.StringIO()
    with redirect_stdout(buf):
        func(*args, **kwargs)
    return buf.getvalue()


class TestJsonPayload(unittest.TestCase):
    def testCountsAndJobListExcludeSuccess(self):
        collected = {"jobs": [], "rollup": {}, "hints": [], "jobs_not_read": 0}
        payload = ccf.json_payload(
            SAMPLE_INFO, SAMPLE_STATE, collected, max_failures=10
        )
        self.assertEqual(
            payload["counts"],
            {
                "jobs": 3,
                "failed": 1,
                "blocked": 0,
                "pending": 1,
                "pending_blocked": 0,
                "passed": 1,
            },
        )
        self.assertEqual([j["context"] for j in payload["jobs"]], ["A", "C"])
        self.assertFalse(payload["complete"])


class TestPrintDigest(unittest.TestCase):
    def testReportsCountsAndIncompleteWarning(self):
        text = captured(ccf.print_digest, SAMPLE_INFO, SAMPLE_STATE)
        self.assertIn("PR #42 br @ abc123d", text)
        self.assertIn("1 failed, 1 pending, 1 passed / 3", text)
        self.assertIn("EVENT INCOMPLETE", text)
        self.assertIn("A  u1", text)


class TestPrintJobHeader(unittest.TestCase):
    def testKnownFieldsPrintedUnknownFieldsAreQuestionMarks(self):
        env = {
            "CIVET_RECIPE_NAME": "My Recipe",
            "CIVET_PR_NUM": "99",
            "CIVET_EVENT_CAUSE": "Pull request",
        }
        text = captured(ccf.print_job_header, "http://job/2", env)
        self.assertIn("My Recipe  http://job/2", text)
        self.assertIn("pull request 99 (Pull request)", text)
        self.assertIn("head ?:? @ ?", text)
        self.assertIn("base ?:? @ ?", text)

    def testPushEventReportsEventCauseInPlaceOfPullRequest(self):
        # A push to devel carries no CIVET_PR_NUM, so the event cause is
        # all that names the event and the head commit is what to follow it with
        env = {
            "CIVET_RECIPE_NAME": "My Recipe",
            "CIVET_EVENT_CAUSE": "Push devel",
            "CIVET_HEAD_REPO": "libMesh/libmesh",
            "CIVET_HEAD_REF": "devel",
            "CIVET_HEAD_SHA": "0123456789abcdef",
        }
        text = captured(ccf.print_job_header, "http://job/3", env)
        self.assertIn("Push devel", text)
        self.assertNotIn("pull request", text)
        self.assertIn("head libMesh/libmesh:devel @ 0123456789ab", text)


class TestPrintErrorsAndSteps(unittest.TestCase):
    def setUp(self):
        self.steps = [
            ("01_Build", "completed with return code 1\nerror: something broke\n")
        ]
        self.env = {"CIVET_RECIPE_NAME": "My Recipe"}

    def testPrintErrorsShowsDiagnostics(self):
        collected = ccf.collect_one_job("http://job/2", self.steps, self.env, 5, 5)
        text = captured(ccf.print_errors, collected, max_failures=10, label_jobs=False)
        self.assertIn("01_Build (exit 1, build)", text)
        self.assertIn("error: something broke", text)

    def testPrintStepListShowsSizeInKB(self):
        text = captured(ccf.print_step_list, "http://job/2", self.steps)
        self.assertIn("steps in http://job/2", text)
        self.assertIn("01_Build", text)

    def testPrintStepLogTailWhenNoPattern(self):
        text = captured(ccf.print_step_log, self.steps, "01_Build", 100, None, 0)
        self.assertIn("last 2 of 2 lines", text)
        self.assertIn("error: something broke", text)

    def testPrintStepLogGrepsAndMarksGaps(self):
        text = captured(ccf.print_step_log, self.steps, "01_Build", 100, "broke", 0)
        self.assertIn("1 of 2 lines match 'broke'", text)
        self.assertIn("error: something broke", text)

    def testPrintStepLogUnknownStepRaises(self):
        with self.assertRaises(ccf.GitHubError):
            ccf.print_step_log(self.steps, "no-such-step", 100, None, 0)


def completed(returncode, stdout="", stderr=""):
    return subprocess.CompletedProcess(
        args=[], returncode=returncode, stdout=stdout, stderr=stderr
    )


class TestGh(unittest.TestCase):
    def testSuccessReturnsStdout(self):
        with mock.patch("subprocess.run", return_value=completed(0, "hello\n")):
            self.assertEqual(ccf.gh(["foo"]), "hello\n")

    def testNonZeroExitRaisesGitHubError(self):
        with mock.patch("subprocess.run", return_value=completed(1, "", "boom")):
            with self.assertRaisesRegex(ccf.GitHubError, "boom"):
                ccf.gh(["foo"])

    def testMissingGhRaisesGitHubError(self):
        with mock.patch("subprocess.run", side_effect=FileNotFoundError()):
            with self.assertRaisesRegex(ccf.GitHubError, "not found on PATH"):
                ccf.gh(["foo"])

    def testGhJsonParsesOrReturnsNoneWhenEmpty(self):
        with mock.patch.object(ccf, "gh", return_value='{"a": 1}\n'):
            self.assertEqual(ccf.gh_json(["x"]), {"a": 1})
        with mock.patch.object(ccf, "gh", return_value=""):
            self.assertIsNone(ccf.gh_json(["x"]))


class TestCurrentBranchAndRepoFromRemote(unittest.TestCase):
    def testCurrentBranch(self):
        with mock.patch("subprocess.run", return_value=completed(0, "mybranch\n")):
            self.assertEqual(ccf.current_branch(), "mybranch")

    def testCurrentBranchNotAGitRepoRaises(self):
        with mock.patch("subprocess.run", return_value=completed(1)):
            with self.assertRaises(ccf.GitHubError):
                ccf.current_branch()

    def testRepoFromRemoteParsesSlug(self):
        with mock.patch(
            "subprocess.run",
            return_value=completed(0, "git@github.com:libMesh/libmesh.git\n"),
        ):
            self.assertEqual(ccf.repo_from_remote("origin"), "libMesh/libmesh")

    def testRepoFromRemoteMissingRemoteReturnsNone(self):
        with mock.patch("subprocess.run", return_value=completed(1)):
            self.assertIsNone(ccf.repo_from_remote("nope"))


class TestResolveRepo(unittest.TestCase):
    def testExplicitRepoWins(self):
        self.assertEqual(ccf.resolve_repo("owner/name", None), "owner/name")

    def testNamedRemoteUsedBeforeCandidates(self):
        with mock.patch.object(
            ccf, "repo_from_remote", return_value="owner/fromremote"
        ):
            self.assertEqual(ccf.resolve_repo(None, "up"), "owner/fromremote")

    def testCandidatesTriedInOrderUntilOneHits(self):
        seen = []

        def fake(remote):
            seen.append(remote)
            return "owner/x" if remote == "upstream" else None

        with mock.patch.object(ccf, "repo_from_remote", side_effect=fake):
            self.assertEqual(ccf.resolve_repo(None, None), "owner/x")
        self.assertEqual(seen, ["up", "upstream"])

    def testNoneFoundReturnsNone(self):
        with mock.patch.object(ccf, "repo_from_remote", return_value=None):
            self.assertIsNone(ccf.resolve_repo(None, None))


class TestResolvePrAndTarget(unittest.TestCase):
    def testExplicitPrResolved(self):
        with mock.patch.object(ccf, "gh_json", return_value={"number": 5}):
            self.assertEqual(ccf.resolve_pr("o/n", 5), {"number": 5})

    def testExplicitPrNotFoundRaises(self):
        with mock.patch.object(ccf, "gh_json", return_value=None):
            with self.assertRaises(ccf.GitHubError):
                ccf.resolve_pr("o/n", 5)

    def testDefaultsToBranchsOpenPr(self):
        with (
            mock.patch.object(ccf, "current_branch", return_value="mybr"),
            mock.patch.object(ccf, "gh_json", return_value=[{"number": 7}]),
        ):
            self.assertEqual(ccf.resolve_pr("o/n", None), {"number": 7})

    def testNoOpenPrForBranchRaises(self):
        with (
            mock.patch.object(ccf, "current_branch", return_value="mybr"),
            mock.patch.object(ccf, "gh_json", return_value=[]),
        ):
            with self.assertRaises(ccf.GitHubError):
                ccf.resolve_pr("o/n", None)

    def testResolveTargetShaSkipsPrLookup(self):
        info = ccf.resolve_target("o/n", None, "deadbeef")
        self.assertEqual(info["headRefOid"], "deadbeef")
        self.assertIsNone(info["number"])

    def testResolveTargetPrDelegatesToResolvePr(self):
        with mock.patch.object(ccf, "resolve_pr", return_value={"number": 1}):
            self.assertEqual(ccf.resolve_target("o/n", 1, None), {"number": 1})


class FakeUrlResponse:
    """Minimal context-manager stand-in for urllib.request.urlopen's result."""

    def __init__(self, data):
        self.data = data

    def read(self):
        return self.data

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class TestFetchState(unittest.TestCase):
    def testCombinesRollupStateAndPaginatedStatuses(self):
        status_line = json.dumps({"context": "A", "state": "failure"}) + "\n"
        with (
            mock.patch.object(ccf, "gh_json", return_value={"state": "failure"}),
            mock.patch.object(ccf, "gh", return_value=status_line),
        ):
            state = ccf.fetch_state("o/n", "sha")
        self.assertEqual(state["state"], "failure")
        self.assertEqual(state["statuses"], [{"context": "A", "state": "failure"}])


class TestFetchJobSteps(unittest.TestCase):
    def testExtractsEachTarMemberAsAStep(self):
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tar:
            data = b"hello world\n"
            info = tarfile.TarInfo(name="job/01_Build")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))

        with mock.patch(
            "urllib.request.urlopen", return_value=FakeUrlResponse(buf.getvalue())
        ):
            steps = ccf.fetch_job_steps("https://civet.inl.gov/job/123/")
        self.assertEqual(steps, [("01_Build", "hello world\n")])


class TestParseArgs(unittest.TestCase):
    def testDefaults(self):
        args = ccf.parse_args([])
        self.assertIsNone(args.pr)
        self.assertIsNone(args.sha)
        self.assertIsNone(args.job)
        self.assertEqual(args.max_failures, ccf.DEFAULT_MAX_FAILURES)
        self.assertFalse(args.as_json)

    def testPrAndShaAreMutuallyExclusive(self):
        with mock.patch("sys.stderr", new_callable=io.StringIO):
            with self.assertRaises(SystemExit):
                ccf.parse_args(["--pr", "5", "--sha", "abc"])


JOB_STEPS = [("01_Build", "completed with return code 1\nerror: something broke\n")]


class TestMainJobDispatch(unittest.TestCase):
    def testListSteps(self):
        with mock.patch.object(ccf, "fetch_job_steps", return_value=JOB_STEPS):
            text = captured(ccf.main, ["--job", "http://job/2", "--list-steps"])
        self.assertIn("steps in http://job/2", text)

    def testStep(self):
        with mock.patch.object(ccf, "fetch_job_steps", return_value=JOB_STEPS):
            text = captured(ccf.main, ["--job", "http://job/2", "--step", "01_Build"])
        self.assertIn("error: something broke", text)

    def testPlainJobReportsHeaderAndErrors(self):
        with mock.patch.object(ccf, "fetch_job_steps", return_value=JOB_STEPS):
            text = captured(ccf.main, ["--job", "http://job/2"])
        self.assertIn("http://job/2", text)
        self.assertIn("error: something broke", text)

    def testUnreadableJobRaisesGitHubError(self):
        with mock.patch.object(ccf, "fetch_job_steps", side_effect=RuntimeError("403")):
            with self.assertRaises(ccf.GitHubError):
                ccf.main(["--job", "http://job/2"])


class TestMainPrDispatch(unittest.TestCase):
    def setUp(self):
        self.info = {"number": 5, "headRefName": "br", "headRefOid": "deadbeef123"}
        self.state = {"state": "success", "statuses": []}

    def testDigestOnlyWhenNoErrorsOrJson(self):
        with (
            mock.patch.object(ccf, "resolve_repo", return_value="o/n"),
            mock.patch.object(ccf, "resolve_target", return_value=self.info),
            mock.patch.object(ccf, "fetch_state", return_value=self.state),
            mock.patch.object(ccf, "collect_errors") as collect,
        ):
            text = captured(ccf.main, ["--pr", "5"])
        collect.assert_not_called()
        self.assertIn("PR #5 br @ deadbee", text)

    def testJsonImpliesErrorsAndEmitsValidJson(self):
        collected = {"jobs": [], "rollup": {}, "hints": [], "jobs_not_read": 0}
        with (
            mock.patch.object(ccf, "resolve_repo", return_value="o/n"),
            mock.patch.object(ccf, "resolve_target", return_value=self.info),
            mock.patch.object(ccf, "fetch_state", return_value=self.state),
            mock.patch.object(ccf, "collect_errors", return_value=collected) as collect,
        ):
            text = captured(ccf.main, ["--pr", "5", "--json"])
        collect.assert_called_once()
        self.assertEqual(json.loads(text)["pr"], 5)

    def testNoRepoFoundRaisesGitHubError(self):
        with mock.patch.object(ccf, "resolve_repo", return_value=None):
            with self.assertRaises(ccf.GitHubError):
                ccf.main([])


if __name__ == "__main__":
    unittest.main()
