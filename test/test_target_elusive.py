#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "8"
import numpy as np
import polars as pl
from polars.testing import assert_frame_equal
from bird_tool_utils import in_tempdir
from binchicken.workflow.scripts.target_elusive import get_clusters, pipeline, streaming_pipeline

SAMPLE_DISTANCES_COLUMNS = {
    "query_name": str,
    "match_name": str,
    "jaccard": float,
}

CLUSTERS_COLUMNS = {
    "samples": str,
}

APPRAISE_COLUMNS = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    "found_in": str,
}

TARGETS_COLUMNS = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    "target": str,
}

EDGES_COLUMNS = {
    "style": str,
    "cluster_size": int,
    "samples": str,
    "target_ids": str,
}

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtypes=False, check_row_order=False)

    def test_get_clusters(self):
        sample_distances = pl.DataFrame([
            ["sample_1", "sample_2", 1-0.5],
            ["sample_1", "sample_3", 1-1],
            ["sample_2", "sample_3", 1-0.9],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["sample_1", "sample_2", "sample_3"])

        expected_clusters = pl.DataFrame([
            ["sample_1,sample_2"],
            ["sample_2,sample_3"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(sample_distances, samples)
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_suffix(self):
        sample_distances = pl.DataFrame([
            ["sample_1.1", "sample_2.1", 1-0.5],
            ["sample_1.1", "sample_3.1", 1-1],
            ["sample_2.1", "sample_3.1", 1-0.9],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["sample_1", "sample_2", "sample_3"])

        expected_clusters = pl.DataFrame([
            ["sample_1.1,sample_2.1"],
            ["sample_2.1,sample_3.1"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(sample_distances, samples)
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_suffix_underscore(self):
        sample_distances = pl.DataFrame([
            ["sample_1_1", "sample_2_1", 1-0.5],
            ["sample_1_1", "sample_3_1", 1-1],
            ["sample_2_1", "sample_3_1", 1-0.9],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["sample_1_1", "sample_2_1", "sample_3_1"])

        expected_clusters = pl.DataFrame([
            ["sample_1_1,sample_2_1"],
            ["sample_2_1,sample_3_1"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(sample_distances, samples)
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_suffix_R(self):
        sample_distances = pl.DataFrame([
            ["sample_1_R1", "sample_2_R1", 1-0.5],
            ["sample_1_R1", "sample_3_R1", 1-1],
            ["sample_2_R1", "sample_3_R1", 1-0.9],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["sample_1_R1", "sample_2_R1", "sample_3_R1"])

        expected_clusters = pl.DataFrame([
            ["sample_1_R1,sample_2_R1"],
            ["sample_2_R1,sample_3_R1"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(sample_distances, samples)
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_large_precluster_size(self):
        sample_distances = pl.DataFrame([
            ["sample_1", "sample_2", 1-0.5],
            ["sample_1", "sample_3", 1-1],
            ["sample_2", "sample_3", 1-0.9],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["sample_1", "sample_2", "sample_3"])

        expected_clusters = pl.DataFrame([
            ["sample_1,sample_2"],
            ["sample_2,sample_3"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(sample_distances, samples, PRECLUSTER_SIZE=10)
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_size_three(self):
        sample_distances = pl.DataFrame([
            ["1", "2", 1-0.5],
            ["1", "3", 1-1],
            ["2", "3", 1-0.9],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["1", "2", "3"])

        expected_clusters = pl.DataFrame([
            ["1,2"],
            ["1,3"],
            ["2,3"],
            ["1,2,3"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(
            sample_distances,
            samples,
            PRECLUSTER_SIZE=3,
            MAX_COASSEMBLY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_size_three_of_four(self):
        sample_distances = pl.DataFrame([
            ["1", "2", 1-0.1],
            ["1", "3", 1-0.2],
            ["1", "4", 1-0.4],
            ["2", "3", 1-0.2],
            ["2", "4", 1-0.4],
            ["3", "4", 1-1],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["1", "2", "3", "4"])

        expected_clusters = pl.DataFrame([
            ["1,2"],
            ["1,3"],
            ["1,4"],
            ["2,3"],
            ["2,4"],
            # ["3,4"],
            ["1,2,3"],
            ["1,2,4"],
            # ["1,3,4"],
            # ["2,3,4"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(
            sample_distances,
            samples,
            PRECLUSTER_SIZE=3,
            MAX_COASSEMBLY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_size_three_of_four_missing(self):
        sample_distances = pl.DataFrame([
            ["1", "2", 1-0.1],
            ["1", "3", 1-0.2],
            ["2", "3", 1-0.2],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["1", "2", "3"])

        expected_clusters = pl.DataFrame([
            ["1,2"],
            ["1,3"],
            # ["1,4"],
            ["2,3"],
            # ["2,4"],
            # ["3,4"],
            ["1,2,3"],
            # ["1,2,4"],
            # ["1,3,4"],
            # ["2,3,4"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(
            sample_distances,
            samples,
            PRECLUSTER_SIZE=3,
            MAX_COASSEMBLY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_size_three_of_four_balanced(self):
        sample_distances = pl.DataFrame([
            ["1", "2", 1-0.2],
            ["1", "3", 1-0.3],
            ["1", "4", 1-0.4],
            ["2", "3", 1-1],
            ["2", "4", 1-0.2],
            ["3", "4", 1-0.1],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["1", "2", "3", "4"])

        expected_clusters = pl.DataFrame([
            ["1,2"],
            ["1,3"],
            # ["1,4"],
            # ["2,3"],
            ["2,4"],
            ["3,4"],
            ["1,2,3"],
            ["1,2,4"],
            ["1,3,4"],
            ["2,3,4"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(
            sample_distances,
            samples,
            PRECLUSTER_SIZE=3,
            MAX_COASSEMBLY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_size_three_of_four_balanced_missing(self):
        sample_distances = pl.DataFrame([
            ["1", "2", 1-0.2],
            ["1", "3", 1-0.3],
            ["1", "4", 1-0.4],
            ["2", "4", 1-0.2],
            ["3", "4", 1-0.1],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["1", "2", "3", "4"])

        expected_clusters = pl.DataFrame([
            ["1,2"],
            ["1,3"],
            # ["1,4"],
            # ["2,3"],
            ["2,4"],
            ["3,4"],
            ["1,2,3"],
            ["1,2,4"],
            ["1,3,4"],
            ["2,3,4"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(
            sample_distances,
            samples,
            PRECLUSTER_SIZE=3,
            MAX_COASSEMBLY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_size_four_of_five(self):
        sample_distances = pl.DataFrame([
            ["1", "2", 1-0.2],
            ["1", "3", 1-0.3],
            ["1", "4", 1-0.4],
            ["1", "5", 1-0.5],
            ["2", "3", 1-0.2],
            ["2", "4", 1-0.3],
            ["2", "5", 1-0.4],
            ["3", "4", 1-0.3],
            ["3", "5", 1-0.3],
            ["4", "5", 1-0.1],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["1", "2", "3", "4", "5"])

        expected_clusters = pl.DataFrame([
            ["1,2"],
            ["1,3"],
            ["1,4"],
            # ["1,5"],
            ["2,3"],
            ["2,4"],
            ["2,5"],
            ["3,4"],
            ["3,5"],
            ["4,5"],
            ["1,2,3"],
            ["1,2,4"],
            ["1,3,4"],
            ["2,3,4"],
            # ["1,2,5"],
            # ["1,3,5"],
            ["2,3,5"],
            # ["1,4,5"],
            ["2,4,5"],
            ["3,4,5"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(
            sample_distances,
            samples,
            PRECLUSTER_SIZE=4,
            MAX_COASSEMBLY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_size_four_of_five_missing(self):
        sample_distances = pl.DataFrame([
            ["1", "2", 1-0.2],
            ["1", "3", 1-0.3],
            ["1", "4", 1-0.4],
            ["2", "3", 1-0.2],
            ["2", "4", 1-0.3],
            ["2", "5", 1-0.4],
            ["3", "4", 1-0.3],
            ["3", "5", 1-0.3],
            ["4", "5", 1-0.1],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["1", "2", "3", "4", "5"])

        expected_clusters = pl.DataFrame([
            ["1,2"],
            ["1,3"],
            ["1,4"],
            # ["1,5"],
            ["2,3"],
            ["2,4"],
            ["2,5"],
            ["3,4"],
            ["3,5"],
            ["4,5"],
            ["1,2,3"],
            ["1,2,4"],
            ["1,3,4"],
            ["2,3,4"],
            # ["1,2,5"],
            # ["1,3,5"],
            ["2,3,5"],
            # ["1,4,5"],
            ["2,4,5"],
            ["3,4,5"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(
            sample_distances,
            samples,
            PRECLUSTER_SIZE=4,
            MAX_COASSEMBLY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_real_world(self):
        sample_distances = pl.DataFrame([
            ["SRR12149290", "SRR10571243", 1-0.16],
            ["ERR3201415", "ERR3220216", 1-0.16],
            ["ERR1414209", "ERR4804028", 1-0.16],
            ["SRR6979552", "SRR15213103", 1-0.07],
            ["SRR12149290", "SRR12352217", 1-0.16],
            ["SRR6979552", "SRR4831657", 1-0.09],
            ["ERR3201415", "SRR11784293", 1-0.09],
            ["SRR4249921", "SRR5207344", 1-0.15],
            ["SRR6979552", "SRR6980357", 1-0.25],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set([
            "SRR12149290", "SRR10571243",
            "ERR3201415", "ERR3220216",
            "ERR1414209", "ERR4804028",
            "SRR6979552", "SRR15213103",
            "SRR12149290", "SRR12352217",
            "SRR6979552", "SRR4831657",
            "ERR3201415", "SRR11784293",
            "SRR4249921", "SRR5207344",
            "SRR6979552", "SRR6980357",
        ])

        expected_clusters = pl.DataFrame([
            ["SRR10571243,SRR12149290"],
            ["ERR3201415,ERR3220216"],
            ["ERR1414209,ERR4804028"],
            ["SRR15213103,SRR6979552"],
            ["SRR12149290,SRR12352217"],
            ["SRR4831657,SRR6979552"],
            ["ERR3201415,SRR11784293"],
            ["SRR4249921,SRR5207344"],
            ["SRR6979552,SRR6980357"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(sample_distances, samples)
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_empty_inputs(self):
        sample_distances = pl.DataFrame([
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["sample_1", "sample_2"])

        expected_clusters = pl.DataFrame([
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(sample_distances, samples)
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_anchored(self):
        sample_distances = pl.DataFrame([
            ["sample_1", "sample_2", 1-0.5],
            ["sample_1", "sample_3", 1-1],
            ["sample_2", "sample_3", 1-0.9],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["sample_1", "sample_2", "sample_3"])
        anchor_samples = set(["sample_1"])

        expected_clusters = pl.DataFrame([
            ["sample_1,sample_2"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(sample_distances, samples, anchor_samples=anchor_samples)
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_anchored_size_four_of_five(self):
        sample_distances = pl.DataFrame([
            ["1", "2", 1-0.2],
            ["1", "3", 1-0.3],
            ["1", "4", 1-0.4],
            ["1", "5", 1-0.5],
            ["2", "3", 1-0.2],
            ["2", "4", 1-0.3],
            ["2", "5", 1-0.4],
            ["3", "4", 1-0.3],
            ["3", "5", 1-0.3],
            ["4", "5", 1-0.1],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set(["1", "2", "3", "4", "5"])
        anchor_samples = set(["1"])

        expected_clusters = pl.DataFrame([
            ["1,2"],
            ["1,3"],
            ["1,4"],
            # ["1,5"],
            # ["2,3"],
            # ["2,4"],
            # ["2,5"],
            # ["3,4"],
            # ["3,5"],
            # ["4,5"],
            ["1,2,3"],
            ["1,2,4"],
            ["1,3,4"],
            # ["2,3,4"],
            # ["1,2,5"],
            # ["1,3,5"],
            # ["2,3,5"],
            # ["1,4,5"],
            # ["2,4,5"],
            # ["3,4,5"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(
            sample_distances,
            samples,
            anchor_samples=anchor_samples,
            PRECLUSTER_SIZE=4,
            MAX_COASSEMBLY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_get_clusters_anchored_real_world(self):
        sample_distances = pl.DataFrame([
            ["SRR12149290", "SRR10571243", 1-0.16],
            ["ERR3201415", "ERR3220216", 1-0.16],
            ["ERR1414209", "ERR4804028", 1-0.16],
            ["SRR6979552", "SRR15213103", 1-0.07],
            ["SRR12149290", "SRR12352217", 1-0.16],
            ["SRR6979552", "SRR4831657", 1-0.09],
            ["ERR3201415", "SRR11784293", 1-0.09],
            ["SRR4249921", "SRR5207344", 1-0.15],
            ["SRR6979552", "SRR6980357", 1-0.25],
        ], orient="row", schema=SAMPLE_DISTANCES_COLUMNS)
        samples = set([
            "SRR12149290", "SRR10571243",
            "ERR3201415", "ERR3220216",
            "ERR1414209", "ERR4804028",
            "SRR6979552", "SRR15213103",
            "SRR12149290", "SRR12352217",
            "SRR6979552", "SRR4831657",
            "ERR3201415", "SRR11784293",
            "SRR4249921", "SRR5207344",
            "SRR6979552", "SRR6980357",
        ])
        anchor_samples = set(["SRR10571243", "SRR11784293", "SRR6979552"])

        expected_clusters = pl.DataFrame([
            ["SRR10571243,SRR12149290"],
            # ["ERR3201415,ERR3220216"],
            # ["ERR1414209,ERR4804028"],
            ["SRR15213103,SRR6979552"],
            # ["SRR12149290,SRR12352217"],
            # ["SRR4831657,SRR6979552"],
            ["ERR3201415,SRR11784293"],
            # ["SRR4249921,SRR5207344"],
            # ["SRR6979552,SRR6980357"],
        ], orient="row", schema=CLUSTERS_COLUMNS)

        observed_clusters = get_clusters(
            sample_distances,
            samples,
            anchor_samples=anchor_samples,
            )
        self.assertDataFrameEqual(expected_clusters, observed_clusters)

    def test_target_elusive(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_suffix(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2.1", "AAA", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1.1", "sample_2.1"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2.1", "AAA", 5, 10, "Root", "0"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1.1,sample_2.1", "0"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_suffix_underscore(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2_1", "AAA", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1_1", "sample_2_1"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1_1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2_1", "AAA", 5, 10, "Root", "0"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1_1,sample_2_1", "0"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_suffix_R(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1_R1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2_R1", "AAA", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1_R1", "sample_2_R1"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1_R1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2_R1", "AAA", 5, 10, "Root", "0"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1_R1,sample_2_R1", "0"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_empty_input(self):
        unbinned = pl.DataFrame([
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2"])

        expected_targets = pl.DataFrame([
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_no_matching_samples(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2.1", "AAA", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_3", "sample_4"])

        expected_targets = pl.DataFrame([
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_low_coverage(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 5, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 4, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 5, "Root", "0"],
            ["S3.1", "sample_2", "AAA", 5, 4, "Root", "0"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_two_targets(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", "1"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0,1"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_multiple_genes(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.2", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.2", "sample_2", "AAB", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.2", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.2", "sample_2", "AAB", 5, 10, "Root", "1"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0,1"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_multiple_genes_same_sequence(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.2", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.2", "sample_2", "AAA", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.2", "sample_1", "AAA", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.2", "sample_2", "AAA", 5, 10, "Root", "1"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0,1"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_three_sample_mix(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_3", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_3", "AAC", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2", "sample_3"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_3", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_3", "AAC", 5, 10, "Root", "2"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0,1"],
            ["match", 2, "sample_1,sample_3", "0"],
            ["match", 2, "sample_2,sample_3", "0"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_skip_samples(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_3", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_3", "AAC", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", "1"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0,1"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_three_way_targets(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 3, 6, "Root", ""],
            ["S3.1", "sample_1", "AAB", 2, 4, "Root", ""],
            ["S3.1", "sample_2", "AAA", 3, 6, "Root", ""],
            ["S3.1", "sample_2", "AAB", 2, 4, "Root", ""],
            ["S3.1", "sample_3", "AAA", 2, 4, "Root", ""],
            ["S3.1", "sample_3", "AAB", 2, 4, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2", "sample_3"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 3, 6, "Root", "0"],
            ["S3.1", "sample_1", "AAB", 2, 4, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 3, 6, "Root", "0"],
            ["S3.1", "sample_2", "AAB", 2, 4, "Root", "1"],
            ["S3.1", "sample_3", "AAA", 2, 4, "Root", "0"],
            ["S3.1", "sample_3", "AAB", 2, 4, "Root", "1"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0"],
            ["pool", 3, "sample_1,sample_2,sample_3", "0,1"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples, MAX_COASSEMBLY_SAMPLES=3)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_four_way_targets(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 2, 4, "Root", ""],
            ["S3.1", "sample_1", "AAC", 1, 3, "Root", ""],
            ["S3.1", "sample_2", "AAA", 2, 4, "Root", ""],
            ["S3.1", "sample_2", "AAB", 2, 4, "Root", ""],
            ["S3.1", "sample_2", "AAC", 1, 3, "Root", ""],
            ["S3.1", "sample_3", "AAA", 2, 4, "Root", ""],
            ["S3.1", "sample_3", "AAB", 2, 4, "Root", ""],
            ["S3.1", "sample_3", "AAC", 1, 3, "Root", ""],
            ["S3.1", "sample_4", "AAB", 2, 4, "Root", ""],
            ["S3.1", "sample_4", "AAC", 1, 3, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2", "sample_3", "sample_4"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 2, 4, "Root", "0"],
            ["S3.1", "sample_1", "AAC", 1, 3, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 2, 4, "Root", "0"],
            ["S3.1", "sample_2", "AAB", 2, 4, "Root", "2"],
            ["S3.1", "sample_2", "AAC", 1, 3, "Root", "1"],
            ["S3.1", "sample_3", "AAA", 2, 4, "Root", "0"],
            ["S3.1", "sample_3", "AAB", 2, 4, "Root", "2"],
            ["S3.1", "sample_3", "AAC", 1, 3, "Root", "1"],
            ["S3.1", "sample_4", "AAB", 2, 4, "Root", "2"],
            ["S3.1", "sample_4", "AAC", 1, 3, "Root", "1"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["pool", 3, "sample_1,sample_2,sample_3", "0"],
            ["pool", 3, "sample_2,sample_3,sample_4", "2"],
            ["pool", 4, "sample_1,sample_2,sample_3,sample_4", "1"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples, MAX_COASSEMBLY_SAMPLES=4)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    @unittest.skip("Benchmarking")
    def test_target_elusive_benchmarking(self):
        genes = ["S3.1", "S3.2", "S3.3", "S3.4", "S3.5", "S3.6", "S3.7", "S3.8", "S3.9", "S3.10"]
        samples = ["sample_" + str(n) for n in range(1, 51)]
        sequences = ["AAA", "AAB", "AAC", "AAD", "AAE", "AAF", "AAG", "AAH", "AAI", "AAJ"]
        unbinned = pl.DataFrame([
            [g, s, n, 2, 5.1, "Root", ""] for g in genes for s in samples for n in sequences
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(samples)

        _, _ = pipeline(unbinned, samples, MAX_COASSEMBLY_SAMPLES=4)
        _, _ = pipeline(unbinned, samples, MAX_COASSEMBLY_SAMPLES=4)
        _, _ = pipeline(unbinned, samples, MAX_COASSEMBLY_SAMPLES=4)
        _, _ = pipeline(unbinned, samples, MAX_COASSEMBLY_SAMPLES=4)
        _, _ = pipeline(unbinned, samples, MAX_COASSEMBLY_SAMPLES=4)

        # Cross join, recursive time: 247.833s
        # Pooling time: 1.979s

    def test_target_elusive_target_taxa(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", "0"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", "0"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples, TAXA_OF_INTEREST="p__Planctomycetota")
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_target_taxa_multiple(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root; d__Bacteria; p__Pseudomonadota", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", "0"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root; d__Bacteria; p__Pseudomonadota", "0"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples, TAXA_OF_INTEREST="p__Planctomycetota|p__Pseudomonadota")
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_single_assembly(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_3", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_3", "AAC", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)
        samples = set(["sample_1", "sample_2", "sample_3"])

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_3", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_3", "AAC", 5, 10, "Root", "2"],
        ], orient="row", schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0,1"],
            ["match", 2, "sample_1,sample_3", "0"],
            ["match", 2, "sample_2,sample_3", "0"],
        ], orient="row", schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, samples, MAX_COASSEMBLY_SAMPLES=1)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_preclustered_four_way(self):
        with in_tempdir():
            unbinned = pl.LazyFrame([
                ["S3.1", "sample_1", "AAA", 2, 4, "Root", ""],
                ["S3.1", "sample_1", "AAC", 1, 3.5, "Root", ""],
                ["S3.1", "sample_2", "AAA", 2, 4, "Root", ""],
                ["S3.1", "sample_2", "AAB", 2, 4, "Root", ""],
                ["S3.1", "sample_2", "AAC", 1, 3.5, "Root", ""],
                ["S3.1", "sample_3", "AAA", 2, 4, "Root", ""],
                ["S3.1", "sample_3", "AAB", 2, 4, "Root", ""],
                ["S3.1", "sample_3", "AAC", 1, 3.5, "Root", ""],
                ["S3.1", "sample_4", "AAB", 2, 4, "Root", ""],
                ["S3.1", "sample_4", "AAC", 1, 3.5, "Root", ""],
            ], orient="row", schema=APPRAISE_COLUMNS)
            samples = set(["sample_1", "sample_2", "sample_3", "sample_4"])
            preclusters = pl.DataFrame([
                ["sample_1,sample_3,sample_4"],
                ["sample_1,sample_2,sample_4"],
                ["sample_1,sample_2,sample_3,sample_4"],
            ], orient="row", schema=CLUSTERS_COLUMNS)

            expected_targets = pl.DataFrame([
                ["S3.1", "sample_1", "AAA", 2, 4, "Root", "3742656936210706615"],
                ["S3.1", "sample_1", "AAC", 1, 3.5, "Root", "4119220645959756749"],
                ["S3.1", "sample_2", "AAA", 2, 4, "Root", "3742656936210706615"],
                ["S3.1", "sample_2", "AAB", 2, 4, "Root", "3678881348843394362"],
                ["S3.1", "sample_2", "AAC", 1, 3.5, "Root", "4119220645959756749"],
                ["S3.1", "sample_3", "AAA", 2, 4, "Root", "3742656936210706615"],
                ["S3.1", "sample_3", "AAB", 2, 4, "Root", "3678881348843394362"],
                ["S3.1", "sample_3", "AAC", 1, 3.5, "Root", "4119220645959756749"],
                ["S3.1", "sample_4", "AAB", 2, 4, "Root", "3678881348843394362"],
                ["S3.1", "sample_4", "AAC", 1, 3.5, "Root", "4119220645959756749"],
            ], orient="row", schema=TARGETS_COLUMNS)
            expected_edges = pl.DataFrame([
                ["match", 3, "sample_1,sample_3,sample_4", "4119220645959756749"],
                ["match", 3, "sample_1,sample_2,sample_4", "4119220645959756749"],
                ["match", 4, "sample_1,sample_2,sample_3,sample_4", "4119220645959756749"],
            ], orient="row", schema=EDGES_COLUMNS)

            streaming_pipeline(
                unbinned,
                samples,
                sample_preclusters=preclusters,
                targets_path="targets.tsv",
                edges_path="edges.tsv",
                MAX_COASSEMBLY_SAMPLES=2,
                CHUNK_SIZE=2,
                )
            observed_targets = pl.read_csv("targets.tsv", schema_overrides=TARGETS_COLUMNS, separator="\t")
            observed_edges = pl.read_csv("edges.tsv", schema_overrides=EDGES_COLUMNS, separator="\t")
            self.assertDataFrameEqual(expected_targets, observed_targets)
            self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_preclustered_single_assembly(self):
        with in_tempdir():
            unbinned = pl.LazyFrame([
                ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
                ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
                ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
                ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
                ["S3.1", "sample_3", "AAA", 5, 10, "Root", ""],
                ["S3.1", "sample_3", "AAC", 5, 10, "Root", ""],
            ], orient="row", schema=APPRAISE_COLUMNS)
            samples = set(["sample_1", "sample_2", "sample_3"])
            preclusters = pl.DataFrame([
                # ["sample_1,sample_2"],
                ["sample_1,sample_3"],
                ["sample_2,sample_3"],
            ], orient="row", schema=CLUSTERS_COLUMNS)

            expected_targets = pl.DataFrame([
                ["S3.1", "sample_1", "AAA", 5, 10, "Root", "3742656936210706615"],
                ["S3.1", "sample_1", "AAB", 5, 10, "Root", "3678881348843394362"],
                ["S3.1", "sample_2", "AAA", 5, 10, "Root", "3742656936210706615"],
                ["S3.1", "sample_2", "AAB", 5, 10, "Root", "3678881348843394362"],
                ["S3.1", "sample_3", "AAA", 5, 10, "Root", "3742656936210706615"],
                ["S3.1", "sample_3", "AAC", 5, 10, "Root", "4119220645959756749"],
            ], orient="row", schema=TARGETS_COLUMNS)
            expected_edges = pl.DataFrame([
                # ["match", 2, "sample_1,sample_2", "3742656936210706615,3678881348843394362"],
                ["match", 2, "sample_1,sample_3", "3742656936210706615"],
                ["match", 2, "sample_2,sample_3", "3742656936210706615"],
            ], orient="row", schema=EDGES_COLUMNS)

            streaming_pipeline(
                unbinned,
                samples,
                sample_preclusters=preclusters,
                targets_path="targets.tsv",
                edges_path="edges.tsv",
                MAX_COASSEMBLY_SAMPLES=1,
                CHUNK_SIZE=2,
                )
            observed_targets = pl.read_csv("targets.tsv", schema_overrides=TARGETS_COLUMNS, separator="\t")
            observed_edges = pl.read_csv("edges.tsv", schema_overrides=EDGES_COLUMNS, separator="\t")
            self.assertDataFrameEqual(expected_targets, observed_targets)
            self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_preclustered_coalesce_target(self):
        with in_tempdir():
            unbinned = pl.LazyFrame([
                ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
                ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
                ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
                ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
                ["S3.1", "sample_3", "AAA", 5, 10, "Root", ""],
                ["S3.1", "sample_3", "AAC", 5, 10, "Root", ""],
            ], orient="row", schema=APPRAISE_COLUMNS)
            samples = set(["sample_1", "sample_2", "sample_3"])
            preclusters = pl.DataFrame([
                ["sample_1,sample_2"],
                ["sample_1,sample_3"],
                ["sample_2,sample_3"],
            ], orient="row", schema=CLUSTERS_COLUMNS)

            expected_targets = pl.DataFrame([
                ["S3.1", "sample_1", "AAA", 5, 10, "Root", "3742656936210706615"],
                ["S3.1", "sample_1", "AAB", 5, 10, "Root", "3678881348843394362"],
                ["S3.1", "sample_2", "AAA", 5, 10, "Root", "3742656936210706615"],
                ["S3.1", "sample_2", "AAB", 5, 10, "Root", "3678881348843394362"],
                ["S3.1", "sample_3", "AAA", 5, 10, "Root", "3742656936210706615"],
                ["S3.1", "sample_3", "AAC", 5, 10, "Root", "4119220645959756749"],
            ], orient="row", schema=TARGETS_COLUMNS)
            expected_edges = pl.DataFrame([
                ["match", 2, "sample_1,sample_2", "3678881348843394362,3742656936210706615"],
                ["match", 2, "sample_1,sample_3", "3742656936210706615"],
                ["match", 2, "sample_2,sample_3", "3742656936210706615"],
            ], orient="row", schema=EDGES_COLUMNS)

            streaming_pipeline(
                unbinned,
                samples,
                sample_preclusters=preclusters,
                targets_path="targets.tsv",
                edges_path="edges.tsv",
                MAX_COASSEMBLY_SAMPLES=1,
                CHUNK_SIZE=2,
                )
            observed_targets = pl.read_csv("targets.tsv", schema_overrides=TARGETS_COLUMNS, separator="\t")
            observed_edges = pl.read_csv("edges.tsv", schema_overrides=EDGES_COLUMNS, separator="\t")
            self.assertDataFrameEqual(expected_targets, observed_targets)
            self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_preclustered_no_targets(self):
        with in_tempdir():
            unbinned = pl.LazyFrame([
                ["S3.1", "sample_1", "AAA", 1, 3, "Root", ""],
                ["S3.2", "sample_1", "AAB", 1, 3, "Root", ""],
                ["S3.1", "sample_2", "AAA", 1, 3, "Root", ""],
                ["S3.2", "sample_2", "AAB", 1, 3, "Root", ""],
            ], orient="row", schema=APPRAISE_COLUMNS)
            samples = set(["sample_1", "sample_2"])
            preclusters = pl.DataFrame([
                ["sample_1,sample_2"],
            ], orient="row", schema=CLUSTERS_COLUMNS)

            expected_targets = pl.DataFrame([
                ["S3.1", "sample_1", "AAA", 1, 3, "Root", "3742656936210706615"],
                ["S3.2", "sample_1", "AAB", 1, 3, "Root", "376802398220639393"],
                ["S3.1", "sample_2", "AAA", 1, 3, "Root", "3742656936210706615"],
                ["S3.2", "sample_2", "AAB", 1, 3, "Root", "376802398220639393"],
            ], orient="row", schema=TARGETS_COLUMNS)
            expected_edges = pl.DataFrame([
            ], orient="row", schema=EDGES_COLUMNS)

            streaming_pipeline(
                unbinned,
                samples,
                sample_preclusters=preclusters,
                targets_path="targets.tsv",
                edges_path="edges.tsv",
                CHUNK_SIZE=2,
                )
            observed_targets = pl.read_csv("targets.tsv", schema_overrides=TARGETS_COLUMNS, separator="\t")
            observed_edges = pl.read_csv("edges.tsv", schema_overrides=EDGES_COLUMNS, separator="\t")
            self.assertDataFrameEqual(expected_targets, observed_targets)
            self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_preclustered_empty(self):
        with in_tempdir():
            unbinned = pl.LazyFrame([
                ["S3.1", "sample_1", "AAA", 2, 4, "Root", ""],
                ["S3.1", "sample_1", "AAC", 1, 3.5, "Root", ""],
                ["S3.1", "sample_2", "AAA", 2, 4, "Root", ""],
                ["S3.1", "sample_2", "AAB", 2, 4, "Root", ""],
                ["S3.1", "sample_2", "AAC", 1, 3.5, "Root", ""],
                ["S3.1", "sample_3", "AAA", 2, 4, "Root", ""],
                ["S3.1", "sample_3", "AAB", 2, 4, "Root", ""],
                ["S3.1", "sample_3", "AAC", 1, 3.5, "Root", ""],
                ["S3.1", "sample_4", "AAB", 2, 4, "Root", ""],
                ["S3.1", "sample_4", "AAC", 1, 3.5, "Root", ""],
            ], orient="row", schema=APPRAISE_COLUMNS)
            samples = set(["sample_1", "sample_2", "sample_3", "sample_4"])
            preclusters = pl.DataFrame([
            ], orient="row", schema=CLUSTERS_COLUMNS)

            expected_targets = pl.DataFrame([
                ["S3.1", "sample_1", "AAA", 2, 4, "Root", "3742656936210706615"],
                ["S3.1", "sample_1", "AAC", 1, 3.5, "Root", "4119220645959756749"],
                ["S3.1", "sample_2", "AAA", 2, 4, "Root", "3742656936210706615"],
                ["S3.1", "sample_2", "AAB", 2, 4, "Root", "3678881348843394362"],
                ["S3.1", "sample_2", "AAC", 1, 3.5, "Root", "4119220645959756749"],
                ["S3.1", "sample_3", "AAA", 2, 4, "Root", "3742656936210706615"],
                ["S3.1", "sample_3", "AAB", 2, 4, "Root", "3678881348843394362"],
                ["S3.1", "sample_3", "AAC", 1, 3.5, "Root", "4119220645959756749"],
                ["S3.1", "sample_4", "AAB", 2, 4, "Root", "3678881348843394362"],
                ["S3.1", "sample_4", "AAC", 1, 3.5, "Root", "4119220645959756749"],
            ], orient="row", schema=TARGETS_COLUMNS)
            expected_edges = pl.DataFrame([
            ], orient="row", schema=EDGES_COLUMNS)

            streaming_pipeline(
                unbinned,
                samples,
                sample_preclusters=preclusters,
                targets_path="targets.tsv",
                edges_path="edges.tsv",
                MAX_COASSEMBLY_SAMPLES=2,
                CHUNK_SIZE=2,
                )
            observed_targets = pl.read_csv("targets.tsv", schema_overrides=TARGETS_COLUMNS, separator="\t")
            observed_edges = pl.read_csv("edges.tsv", schema_overrides=EDGES_COLUMNS, separator="\t")
            self.assertDataFrameEqual(expected_targets, observed_targets)
            self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_preclustered_empty_unbinned(self):
        with in_tempdir():
            unbinned = pl.LazyFrame([
            ], orient="row", schema=APPRAISE_COLUMNS)
            samples = set(["sample_1", "sample_2", "sample_3", "sample_4"])
            preclusters = pl.DataFrame([
            ], orient="row", schema=CLUSTERS_COLUMNS)

            expected_targets = pl.DataFrame([
            ], orient="row", schema=TARGETS_COLUMNS)
            expected_edges = pl.DataFrame([
            ], orient="row", schema=EDGES_COLUMNS)

            streaming_pipeline(
                unbinned,
                samples,
                sample_preclusters=preclusters,
                targets_path="targets.tsv",
                edges_path="edges.tsv",
                MAX_COASSEMBLY_SAMPLES=2,
                CHUNK_SIZE=2,
                )
            observed_targets = pl.read_csv("targets.tsv", schema_overrides=TARGETS_COLUMNS, separator="\t")
            observed_edges = pl.read_csv("edges.tsv", schema_overrides=EDGES_COLUMNS, separator="\t")
            self.assertDataFrameEqual(expected_targets, observed_targets)
            self.assertDataFrameEqual(expected_edges, observed_edges)


if __name__ == '__main__':
    unittest.main()
