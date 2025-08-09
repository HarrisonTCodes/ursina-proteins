import unittest
from os import path

from ursina_proteins.protein import Protein, parse_segments


class TestFillSegments(unittest.TestCase):
    def test_basic_usage(self):
        segments = [(5, 10), (14, 15)]
        size = 20
        expected = {
            "a": [(5, 10), (14, 15)],
            "b": [(0, 5), (10, 14), (15, 20)],
        }
        self.assertEqual(parse_segments(segments, size, "a", "b"), expected)

    def test_start_with_segment(self):
        segments = [(0, 10)]
        size = 20
        expected = {"a": [(0, 10)], "b": [(10, 20)]}
        self.assertEqual(parse_segments(segments, size, "a", "b"), expected)

    def test_empty_segments(self):
        segments = []
        size = 12
        expected = {"a": [], "b": [(0, 12)]}
        self.assertEqual(parse_segments(segments, size, "a", "b"), expected)


class TestProteinInit(unittest.TestCase):
    def test_simple_pdb(self):
        base_dir = path.dirname(__file__)
        Protein(path.join(base_dir, "..", "assets", "insulin.pdb"))
