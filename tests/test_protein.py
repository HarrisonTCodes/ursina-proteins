import unittest

from ursina_proteins.protein import fill_segments


class TestFillSegments(unittest.TestCase):

    def test_basic_usage(self):
        segments = [(5, 10), (14, 15)]
        size = 20
        expected = [(0, 4), (5, 10), (11, 13), (14, 15), (16, 20)]
        self.assertEqual(fill_segments(segments, size), expected)

    def test_start_with_segment(self):
        segments = [(0, 10)]
        size = 20
        expected = [(0, 10), (11, 20)]
        self.assertEqual(fill_segments(segments, size), expected)

    def test_empty_segments(self):
        segments = []
        size = 12
        expected = [(0, 12)]
        self.assertEqual(fill_segments(segments, size), expected)
