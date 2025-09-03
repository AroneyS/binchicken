import os
import io
import logging
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import binchicken.binchicken as bc

class DummyProc:
    def __init__(self, stdout_text: str, stderr_text: str = "", returncode: int = 0):
        self._stdout_text = stdout_text
        self._stderr_text = stderr_text
        self.stdout = io.StringIO(stdout_text)
        self.stderr = io.StringIO(stderr_text)
        self.returncode = returncode

    def wait(self):
        return self.returncode

    def poll(self):
        return self.returncode

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if self.stdout:
            self.stdout.close()
        if self.stderr:
            self.stderr.close()


def make_config(path: Path) -> str:
    cfg = path / "config.yaml"
    cfg.write_text("a: 1\n")
    return str(cfg)


class Tests(unittest.TestCase):
    def test_parse_snakemake_errors_dedup_and_order(self):
        text = (
            "Building DAG of jobs...\n"
            "Error in rule aviary_assemble:\n"
            "  jobid: 1\n"
            "Error in rule aviary_recover:\n"
            "  jobid: 2\n"
            "Error in rule aviary_assemble:\n"
            "  jobid: 3\n"
        )
        rules = bc.parse_snakemake_errors(text)
        self.assertEqual(rules, ["aviary_assemble", "aviary_recover"])  # dedup preserves first-seen order

    def test_logs_dir_root_for_workflow(self):
        root = bc.logs_dir_root_for_workflow("/tmp/out", "coassemble.smk")
        self.assertEqual(root, "/tmp/out/coassemble/logs")
        root2 = bc.logs_dir_root_for_workflow("/tmp/out", "evaluate.smk")
        self.assertEqual(root2, "/tmp/out/evaluate/logs")

    def test_find_logs_for_rule_parses_snakefile_and_globs(self):
        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            wid = "UNITTESTWID"
            # Create logs as produced by resources:log_path in coassemble.smk for rule aviary_assemble
            log_dir = tmp_path / "coassemble" / "logs" / "aviary" / "coassembly_0_assemble" / wid
            log_dir.mkdir(parents=True, exist_ok=True)
            (log_dir / "attempt1.log").write_text("first attempt\n")
            (log_dir / "attempt2.log").write_text("second attempt\n")

            found = bc.find_logs_for_rule("aviary_assemble", "coassemble.smk", str(tmp_path), wid)
            # Should find both, with higher attempt first
            self.assertTrue(found)
            self.assertTrue(str(log_dir / "attempt2.log") in found[0])
            self.assertTrue(str(log_dir / "attempt1.log") in found[-1])

            # When logs directory missing, should return empty list
            with TemporaryDirectory() as tmp2:
                found2 = bc.find_logs_for_rule("aviary_assemble", "coassemble.smk", tmp2, wid)
                self.assertEqual(found2, [])

    def test_run_workflow_dumps_logs_and_exits(self):
        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            # Create logs in the structure produced by resources:log_path
            wid = "TESTWID"
            log1 = tmp_path / "coassemble" / "logs" / "aviary" / "coassembly_0_assemble" / wid / "attempt1.log"
            log2 = tmp_path / "coassemble" / "logs" / "aviary" / "coassembly_1_recover" / wid / "attempt1.log"
            log1.parent.mkdir(parents=True, exist_ok=True)
            log2.parent.mkdir(parents=True, exist_ok=True)
            log1.write_text("AAA\nline2\n")
            log2.write_text("BBB\n")

            smk_output = (
                "Building DAG of jobs...\n"
                "Error in rule aviary_assemble:\n"
                "    jobid: 12\n"
                "    output: assembly\n"
                "    shell:\n"
                "        some shell command\n"
                "        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)\n"
                "\n"
                "Error executing rule download_read on cluster (jobid: 592, external: 10161983, jobscript: snakejob.assemble.592.sh). For error details see the cluster log and the log files of the involved rule(s).\n"
                "Trying to restart job 592.\n"
                "\n\n"
                "Error in rule aviary_recover:\n"
                "    jobid: 13\n"
                "    output: recovery\n"
                "    shell:\n"
                "        another shell command\n"
                "        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)\n"
            )

            def fake_popen(cmd, shell, stdout, stderr, text, encoding, errors, bufsize):
                return DummyProc(smk_output, stderr_text="", returncode=1)

            cfg = make_config(tmp_path)

            stdout_buf = io.StringIO()
            stderr_buf = io.StringIO()

            # Capture logging
            log_stream = io.StringIO()
            handler = logging.StreamHandler(log_stream)
            root_logger = logging.getLogger()
            root_level = root_logger.level
            root_logger.setLevel(logging.ERROR)
            root_logger.addHandler(handler)

            try:
                with patch.object(bc.subprocess, "Popen", fake_popen), \
                     patch.object(bc, "workflow_identifier", wid, create=True), \
                     patch.object(bc, "load_configfile", lambda p: None, create=True), \
                     redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
                    with self.assertRaises(SystemExit) as cm:
                        bc.run_workflow(
                            config=cfg,
                            workflow="coassemble.smk",
                            output_dir=str(tmp_path),
                            cores=1,
                            dryrun=False,
                            profile=None,
                            local_cores=1,
                            retries=None,
                            snakemake_args="",
                        )

                self.assertEqual(cm.exception.code, 1)

            finally:
                root_logger.removeHandler(handler)
                root_logger.setLevel(root_level)

            out = stdout_buf.getvalue()
            err = stderr_buf.getvalue()

            # Original snakemake output echoed to stdout
            self.assertIn("Error in rule aviary_assemble", out)
            self.assertIn("Error in rule aviary_recover", out)

            # Log contents dumped to stderr with banners (deduplicated)
            begin1 = f"===== BEGIN LOG ({wid}): {log1}"
            begin2 = f"===== BEGIN LOG ({wid}): {log2}"
            self.assertIn(begin1, err)
            self.assertEqual(err.count(begin1), 1)
            self.assertIn("AAA", err)
            self.assertIn("line2", err)
            self.assertIn(f"===== END LOG ({wid}): {log1}", err)

            self.assertIn(begin2, err)
            self.assertEqual(err.count(begin2), 1)
            self.assertIn("BBB", err)
            self.assertIn(f"===== END LOG ({wid}): {log2}", err)

            # Error logs emitted for each rule/log
            logs = log_stream.getvalue()
            self.assertIn("Rule failed: aviary_assemble; log:", logs)
            self.assertIn("Rule failed: aviary_recover; log:", logs)

    def test_run_workflow_handles_missing_logs(self):
        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            wid = "TESTWID"
            # No logs created under logs dir; parser should report that none were found
            smk_output = "Error in rule build:\n"

            def fake_popen(cmd, shell, stdout, stderr, text, encoding, errors, bufsize):
                return DummyProc(smk_output, stderr_text="", returncode=1)

            cfg = make_config(tmp_path)
            stdout_buf = io.StringIO()
            stderr_buf = io.StringIO()

            # Capture logging
            log_stream = io.StringIO()
            handler = logging.StreamHandler(log_stream)
            root_logger = logging.getLogger()
            root_level = root_logger.level
            root_logger.setLevel(logging.ERROR)
            root_logger.addHandler(handler)

            try:
                with patch.object(bc.subprocess, "Popen", fake_popen), \
                     patch.object(bc, "workflow_identifier", wid, create=True), \
                     patch.object(bc, "load_configfile", lambda p: None, create=True), \
                     redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
                    with self.assertRaises(SystemExit) as cm:
                        bc.run_workflow(
                            config=cfg,
                            workflow="coassemble.smk",
                            output_dir=str(tmp_path),
                            cores=1,
                            dryrun=False,
                            profile=None,
                            local_cores=1,
                            retries=None,
                            snakemake_args="",
                        )
            finally:
                root_logger.removeHandler(handler)
                root_logger.setLevel(root_level)

            self.assertEqual(cm.exception.code, 1)

            out = stdout_buf.getvalue()
            err = stderr_buf.getvalue()
            self.assertIn("Error in rule build", out)
            # Logging should report no logs found under logs dir
            logs = log_stream.getvalue()
            self.assertIn("no log files found under", logs)

    def test_run_workflow_no_parsed_errors_exits(self):
        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            smk_output = "some failure text without the pattern"

            def fake_popen(cmd, shell, stdout, stderr, text, encoding, errors, bufsize):
                return DummyProc(smk_output, stderr_text="", returncode=1)

            cfg = make_config(tmp_path)
            stdout_buf = io.StringIO()
            stderr_buf = io.StringIO()

            with patch.object(bc.subprocess, "Popen", fake_popen), \
                 patch.object(bc, "load_configfile", lambda p: None, create=True), \
                 redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
                with self.assertRaises(SystemExit) as cm:
                    bc.run_workflow(
                        config=cfg,
                        workflow="coassemble.smk",
                        output_dir=str(tmp_path),
                        cores=1,
                        dryrun=False,
                        profile=None,
                        local_cores=1,
                        retries=None,
                        snakemake_args="",
                    )

            self.assertEqual(cm.exception.code, 1)
            out = stdout_buf.getvalue()
            self.assertIn("some failure text", out)


if __name__ == '__main__':
    unittest.main()
