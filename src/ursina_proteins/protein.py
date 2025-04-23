from hashlib import md5

from Bio.PDB import PDBParser
from ursina import Color, Entity, Mesh, color


class Protein:
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

    def __init__(
        self,
        pdb_filepath: str,
        atoms_thickness: float = 0.2,
        chains_thickness: float = 4,
        element_color_map: dict[str, Color] = dict(),
        *args,
        **kwargs,
    ):
        parser = PDBParser()
        structure = parser.get_structure("protein", pdb_filepath)
        structure_centre_of_mass = structure.center_of_mass()

        # Atoms
        atoms_mesh = Mesh(
            mode="point",
            vertices=[atom.coord for atom in structure.get_atoms()],
            colors=[
                element_color_map.get(
                    atom.element,
                    Protein.ELEMENT_COLORS.get(atom.element, color.rgb(1, 0.7, 0.8)),
                )
                for atom in structure.get_atoms()
            ],
            thickness=atoms_thickness,
        )
        self.atoms_entity = Entity(
            model=atoms_mesh, origin=structure_centre_of_mass, *args, **kwargs
        )

        # Chains
        chains_coords = []
        chains_colors = []
        chains_segments = []
        for chain in structure.get_chains():
            segment_start = len(chains_coords)

            # Coords (vertices)
            chain_coords = [
                atom.coord for atom in chain.get_atoms() if atom.get_id() == "CA"
            ]
            chains_coords.extend(chain_coords)

            # Colors
            hash_value = int(md5(chain.get_id().encode("utf-8")).hexdigest(), 16)
            r = (hash_value >> 16) & 0xFF
            g = (hash_value >> 8) & 0xFF
            b = hash_value & 0xFF
            chain_color = color.rgb(r / 255, g / 255, b / 255)
            chains_colors.extend([chain_color for _ in chain_coords])

            # Segments (triangles)
            chains_segments.extend(
                [
                    (i, i + 1)
                    for i in range(segment_start, segment_start + len(chain_coords) - 1)
                ]
            )

        chains_mesh = Mesh(
            mode="line",
            vertices=chains_coords,
            colors=chains_colors,
            triangles=chains_segments,
            thickness=chains_thickness,
        )
        self.chains_entity = Entity(
            model=chains_mesh, origin=structure_centre_of_mass, *args, **kwargs
        )
