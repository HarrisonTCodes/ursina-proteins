from Bio.PDB import PDBParser
from ursina import Entity, Mesh, color


class Protein(Entity):
    ELEMENT_COLORS = {
        "H": color.rgb(1, 1, 1),
        "C": color.rgb(0.2, 0.2, 0.2),
        "N": color.rgb(0, 0, 1),
        "O": color.rgb(1, 0, 0),
        "S": color.rgb(1, 1, 0),
        "P": color.rgb(1, 0.65, 0),
        "Cl": color.rgb(0, 1, 0),
        "Fe": color.rgb(0.7, 0.45, 0.2),
    }

    def __init__(self, pdb_filepath: str, *args, **kwargs):
        parser = PDBParser()
        structure = parser.get_structure("protein", pdb_filepath)

        model = Mesh(
            mode="point",
            vertices=[atom.coord for atom in structure.get_atoms()],
            colors=[
                Protein.ELEMENT_COLORS.get(atom.element, color.rgb(1, 0.7, 0.8))
                for atom in structure.get_atoms()
            ],
            thickness=0.5,
        )

        super().__init__(model=model, *args, **kwargs)
