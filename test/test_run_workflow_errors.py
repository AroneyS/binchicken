import io
import logging
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import binchicken.binchicken as bc


class DummyCompleted:
    def __init__(self, returncode: int, stdout: str):
        self.returncode = returncode
        self.stdout = stdout


def make_config(path: Path) -> str:
    cfg = path / "config.yaml"
    cfg.write_text("a: 1\n")
    return str(cfg)


class Tests(unittest.TestCase):
    def test_run_workflow_dumps_logs_and_exits(self):
        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            log1 = tmp_path / "ruleA.log"
            log2 = tmp_path / "ruleB.log"
            log1.write_text("AAA\nline2\n")
            log2.write_text("BBB\n")

            smk_output = (
                "Building DAG of jobs...\n"
                "Error in rule aviary_assemble:\n"
                "    jobid: 12\n"
                "    output: assembly\n"
                f"    log: {log1} (check log file(s) for error details)\n"
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
                f"    log: {log2} (check log file(s) for error details)\n"
                "    shell:\n"
                "        another shell command\n"
                "        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)\n"
            )

            def fake_run(cmd, shell, stdout, stderr, text, encoding, errors):
                return DummyCompleted(returncode=1, stdout=smk_output)

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
                with patch.object(bc.subprocess, "run", fake_run), \
                     patch.object(bc, "workflow_identifier", "TESTWID", create=True), \
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
                            cluster_retries=None,
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
            begin1 = f"===== BEGIN LOG (TESTWID): {log1}"
            begin2 = f"===== BEGIN LOG (TESTWID): {log2}"
            self.assertIn(begin1, err)
            self.assertEqual(err.count(begin1), 1)
            self.assertIn("AAA", err)
            self.assertIn("line2", err)
            self.assertIn(f"===== END LOG (TESTWID): {log1}", err)

            self.assertIn(begin2, err)
            self.assertEqual(err.count(begin2), 1)
            self.assertIn("BBB", err)
            self.assertIn(f"===== END LOG (TESTWID): {log2}", err)

            # Error logs emitted for each rule/log
            logs = log_stream.getvalue()
            self.assertIn("Rule failed: aviary_assemble; log:", logs)
            self.assertIn("Rule failed: aviary_recover; log:", logs)

    def test_run_workflow_handles_missing_logs(self):
        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            missing_log = tmp_path / "not_there.log"
            smk_output = f"Error in rule build:\n    log: {missing_log}\n"

            def fake_run(cmd, shell, stdout, stderr, text, encoding, errors):
                return DummyCompleted(returncode=1, stdout=smk_output)

            cfg = make_config(tmp_path)
            stdout_buf = io.StringIO()
            stderr_buf = io.StringIO()

            with patch.object(bc.subprocess, "run", fake_run), \
                 patch.object(bc, "workflow_identifier", "TESTWID", create=True), \
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
                        cluster_retries=None,
                        snakemake_args="",
                    )

            self.assertEqual(cm.exception.code, 1)

            out = stdout_buf.getvalue()
            err = stderr_buf.getvalue()
            self.assertIn("Error in rule build", out)
            self.assertIn(f"LOG NOT FOUND (TESTWID): {missing_log}", err)

    def test_run_workflow_no_parsed_errors_exits(self):
        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            smk_output = "some failure text without the pattern"

            def fake_run(cmd, shell, stdout, stderr, text, encoding, errors):
                return DummyCompleted(returncode=1, stdout=smk_output)

            cfg = make_config(tmp_path)
            stdout_buf = io.StringIO()
            stderr_buf = io.StringIO()

            with patch.object(bc.subprocess, "run", fake_run), \
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
                        cluster_retries=None,
                        snakemake_args="",
                    )

            self.assertEqual(cm.exception.code, 1)
            out = stdout_buf.getvalue()
            self.assertIn("some failure text", out)


if __name__ == '__main__':
    unittest.main()
