from Bio.PDB import PDBParser
from ursina import Entity, Mesh


class Protein(Entity):
    def __init__(self, pdb_filepath: str, *args, **kwargs):
        parser = PDBParser()
        structure = parser.get_structure("protein", pdb_filepath)

        atoms = [atom.coord for atom in structure.get_atoms()]
        model = Mesh(mode="point", vertices=atoms)

        super().__init__(model=model, *args, **kwargs)
