import unittest
from os import path

from ursina_proteins.protein import Protein


class TestProteinInit(unittest.TestCase):
    def test_simple_pdb(self):
        base_dir = path.dirname(__file__)
        Protein(path.join(base_dir, "..", "assets", "insulin.pdb"))

    def test_simple_cif(self):
        base_dir = path.dirname(__file__)
        Protein(path.join(base_dir, "..", "assets", "insulin.cif"), legacy_pdb=False)
