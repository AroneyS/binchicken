#!/usr/bin/env python3

import unittest
import os
from bird_tool_utils import in_tempdir
import extern
from snakemake.io import load_configfile

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_conda = os.path.join(path_to_data,'.conda')

METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")
MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
MOCK_COASSEMBLIES = ' '.join([os.path.join(MOCK_COASSEMBLE, "coassembly_0")])

class Tests(unittest.TestCase):
    def test_evaluate(self):
        with in_tempdir():
            cmd = (
                f"cockatoo evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--checkm-version 2 "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            summary_path = os.path.join("test", "evaluate", "evaluate", "summary_stats.tsv")
            self.assertTrue(os.path.exists(summary_path))

            expected = "\n".join(
                [
                    "\t".join([
                        "coassembly",
                        "statistic",
                        "missed",
                        "recovered",
                        "total",
                        "recovered_percent",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "sequences",
                        "1",
                        "1",
                        "2",
                        "50",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "bins",
                        "0",
                        "1",
                        "1",
                        "100",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "taxonomy",
                        "1",
                        "1",
                        "2",
                        "50",
                    ]),
                    ""
                ]
            )
            with open(summary_path) as f:
                self.assertEqual(expected, f.read())

    def test_evaluate_default_config(self):
        with in_tempdir():
            cmd = (
                f"cockatoo evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
            )
            extern.run(cmd)

            config = load_configfile(os.path.join("test", "config.yaml"))
            self.assertEqual(config["max_threads"], 8)
            self.assertEqual(config["checkm_version"], 2)
            self.assertEqual(config["min_completeness"], 70)
            self.assertEqual(config["max_contamination"], 10)


if __name__ == '__main__':
    unittest.main()
