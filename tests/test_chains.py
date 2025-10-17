import unittest

from ursina_proteins.chains import ChainsEntity


class TestParseSegments(unittest.TestCase):
    def test_basic_usage(self):
        segments = [(5, 10), (14, 15)]
        size = 20
        expected = {
            "helices": [(5, 10), (14, 15)],
            "coils": [(0, 5), (10, 14), (15, 20)],
        }
        self.assertEqual(ChainsEntity.parse_segments(segments, size), expected)

    def test_start_with_segment(self):
        segments = [(0, 10)]
        size = 20
        expected = {"helices": [(0, 10)], "coils": [(10, 20)]}
        self.assertEqual(ChainsEntity.parse_segments(segments, size), expected)

    def test_empty_segments(self):
        segments = []
        size = 12
        expected = {"helices": [], "coils": [(0, 12)]}
        self.assertEqual(ChainsEntity.parse_segments(segments, size), expected)
