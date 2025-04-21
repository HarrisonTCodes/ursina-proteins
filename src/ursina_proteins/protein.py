from Bio.PDB import PDBParser
from ursina import Entity


class Protein(Entity):
    def __init__(self, pdb_filepath: str, *args, **kwargs):
        super().__init__(model="cube", *args, **kwargs)
        parser = PDBParser()
        self.structure = parser.get_structure("protein", pdb_filepath)
        print("Center of mass:", self.structure.center_of_mass())
