from Bio.PDB import PDBParser
from ursina import Entity, Mesh, color


class Protein(Entity):
    def __init__(self, pdb_filepath: str, *args, **kwargs):
        parser = PDBParser()
        structure = parser.get_structure("protein", pdb_filepath)

        atoms = [atom.coord for atom in structure.get_atoms()]
        model = Mesh(
            mode="point",
            vertices=atoms,
            colors=[color.random_color() for _ in atoms],
            thickness=0.5,
        )

        super().__init__(model=model, *args, **kwargs)
