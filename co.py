#!/usr/bin/env python3

__author__ = "Samuel Aroney"

import argparse
import logging
import bird_tool_utils as btu

main_parser = btu.BirdArgparser(program="Co",
    examples = {
        "coassemble": [
            btu.Example(
                "coassemble reads and bins",
                "co coassemble --forward reads.1.fq --reverse reads.2.fq --output output_dir"
            ),
            btu.Example(
                "coassemble archive otu tables and bins",
                "co coassemble --singlem-gzip-archives reads.singlem.json.gz --output output_dir"
            )
        ],
        "evaluate": [
            btu.Example(
                "evaluate a completed coassembly",
                "co evaluate --coassembly-output coassembly_dir --output output_dir"
            )
        ]
    }
    )

coassemble_parser = main_parser.new_subparser("coassemble", "Coassemble reads into contigs and bin")

evaluate_parser = main_parser.new_subparser("evaluate", "Evaluate coassembled bins")

args = main_parser.parse_the_args()
