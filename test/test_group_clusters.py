#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from ibis.workflow.scripts.group_clusters import edges_processing, targets_processing, cluster_processing

EDGES_COLUMNS={
    "style": str,
    "cluster_size": int,
    "samples": str,
    "target_ids": str,
    }

TARGETS_COLUMNS={
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    "target": str,
    }

ELUSIVE_CLUSTERS_COLUMNS={
    "samples": str,
    "length": int,
    "total_targets": int,
    "total_size": int,
    "recover_samples": str,
    "coassembly": str,
    }

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtype=False)

    def test_group_edges(self):
        edge1 = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "1,2"],
            ["pool", 3, "sample_1,sample_2,sample_3", "3,4"],
        ], schema=EDGES_COLUMNS)

        edges = [edge1]
        preclusters = ["1"]

        expected = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "1,2"],
            ["pool", 3, "sample_1,sample_2,sample_3", "3,4"],
        ], schema=EDGES_COLUMNS)

        observed = edges_processing(edges, preclusters)
        self.assertDataFrameEqual(expected, observed)

    def test_group_edges_empty_input(self):
        edge1 = pl.DataFrame([
        ], schema=EDGES_COLUMNS)

        edges = [edge1]
        preclusters = ["1"]

        expected = pl.DataFrame([
        ], schema=EDGES_COLUMNS)

        observed = edges_processing(edges, preclusters)
        self.assertDataFrameEqual(expected, observed)

    def test_group_edges_lazy(self):
        edge1 = pl.LazyFrame([
            ["match", 2, "sample_1,sample_2", "1,2"],
            ["pool", 3, "sample_1,sample_2,sample_3", "3,4"],
        ], schema=EDGES_COLUMNS)

        edges = [edge1]
        preclusters = ["1"]

        expected = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "1,2"],
            ["pool", 3, "sample_1,sample_2,sample_3", "3,4"],
        ], schema=EDGES_COLUMNS)

        observed = edges_processing(edges, preclusters).collect()
        self.assertDataFrameEqual(expected, observed)

    def test_group_targets(self):
        target1 = pl.DataFrame([
            ["gene_1", "sample_1", "ATCG", 2, 0.5, "taxonomy_1", "1"],
            ["gene_2", "sample_2", "ATCG", 3, 0.5, "taxonomy_2", "2"],
        ], schema=TARGETS_COLUMNS)

        targets = [target1]
        preclusters = ["1"]

        expected = pl.DataFrame([
            ["gene_1", "sample_1", "ATCG", 2, 0.5, "taxonomy_1", "1"],
            ["gene_2", "sample_2", "ATCG", 3, 0.5, "taxonomy_2", "2"],
        ], schema=TARGETS_COLUMNS)

        observed = targets_processing(targets, preclusters)
        self.assertDataFrameEqual(expected, observed)

    def test_group_targets_empty_input(self):
        target1 = pl.DataFrame([
        ], schema=TARGETS_COLUMNS)

        targets = [target1]
        preclusters = ["1"]

        expected = pl.DataFrame([
        ], schema=TARGETS_COLUMNS)

        observed = targets_processing(targets, preclusters)
        self.assertDataFrameEqual(expected, observed)

    def test_group_targets_lazy(self):
        target1 = pl.LazyFrame([
            ["gene_1", "sample_1", "ATCG", 2, 0.5, "taxonomy_1", "1"],
            ["gene_2", "sample_2", "ATCG", 3, 0.5, "taxonomy_2", "2"],
        ], schema=TARGETS_COLUMNS)

        targets = [target1]
        preclusters = ["1"]

        expected = pl.DataFrame([
            ["gene_1", "sample_1", "ATCG", 2, 0.5, "taxonomy_1", "1"],
            ["gene_2", "sample_2", "ATCG", 3, 0.5, "taxonomy_2", "2"],
        ], schema=TARGETS_COLUMNS)

        observed = targets_processing(targets, preclusters).collect()
        self.assertDataFrameEqual(expected, observed)

    def test_group_clusters(self):
        cluster1 = pl.DataFrame([
            ["sample_1,sample_2", 2, 4, 3, "sample_1,sample_2,sample_3", "coassembly_0"],
            ["sample_2,sample_3", 2, 3, 5, "sample_1,sample_2,sample_3", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)

        clusters = [cluster1]
        preclusters = ["1"]

        expected = pl.DataFrame([
            ["sample_1,sample_2", 2, 4, 3, "sample_1,sample_2,sample_3", "coassembly_0"],
            ["sample_2,sample_3", 2, 3, 5, "sample_1,sample_2,sample_3", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)

        observed = cluster_processing(clusters, preclusters)
        self.assertDataFrameEqual(expected, observed)

    def test_group_clusters_empty_input(self):
        cluster1 = pl.DataFrame([
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)

        clusters = [cluster1]
        preclusters = ["1"]

        expected = pl.DataFrame([
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)

        observed = cluster_processing(clusters, preclusters)
        self.assertDataFrameEqual(expected, observed)

    def test_group_clusters_lazy(self):
        cluster1 = pl.LazyFrame([
            ["sample_1,sample_2", 2, 4, 3, "sample_1,sample_2,sample_3", "coassembly_0"],
            ["sample_2,sample_3", 2, 3, 5, "sample_1,sample_2,sample_3", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)

        clusters = [cluster1]
        preclusters = ["1"]

        expected = pl.DataFrame([
            ["sample_1,sample_2", 2, 4, 3, "sample_1,sample_2,sample_3", "coassembly_0"],
            ["sample_2,sample_3", 2, 3, 5, "sample_1,sample_2,sample_3", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)

        observed = cluster_processing(clusters, preclusters).collect()
        self.assertDataFrameEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
