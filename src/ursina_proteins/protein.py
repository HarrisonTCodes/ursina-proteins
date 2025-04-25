from hashlib import md5

import numpy as np
from Bio.PDB import PDBParser
from scipy.interpolate import make_splrep, splev
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
        chains_smoothness: float = 3,
        *args,
        **kwargs,
    ):
        parser = PDBParser()
        self.structure = parser.get_structure("protein", pdb_filepath)
        structure_center_of_mass = self.structure.center_of_mass()

        self.atoms_entity = Entity(
            model=self.compute_atoms_mesh(atoms_thickness, element_color_map),
            origin=structure_center_of_mass,
            *args,
            **kwargs,
        )

        self.chains_entity = Entity(
            model=self.compute_chains_mesh(chains_thickness, chains_smoothness),
            origin=structure_center_of_mass,
            *args,
            **kwargs,
        )

    def compute_atoms_mesh(
        self, thickness: float, element_color_map: dict[str, Color]
    ) -> Mesh:
        return Mesh(
            mode="point",
            vertices=[atom.coord for atom in self.structure.get_atoms()],
            colors=[
                element_color_map.get(
                    atom.element,
                    Protein.ELEMENT_COLORS.get(atom.element, color.rgb(1, 0.7, 0.8)),
                )
                for atom in self.structure.get_atoms()
            ],
            thickness=thickness,
        )

    def compute_chains_mesh(self, thickness: float, smoothness: float) -> Mesh:
        coords = []
        colors = []
        segments = []

        for chain in self.structure.get_chains():
            segment_start = len(coords)

            # Coords (vertices)
            carbon_alpha_coords = [
                atom.coord for atom in chain.get_atoms() if atom.get_id() == "CA"
            ]

            # Get spline function for each axis
            x, y, z = zip(*carbon_alpha_coords)
            spline_x = make_splrep(range(len(x)), x, s=0)
            spline_y = make_splrep(range(len(y)), y, s=0)
            spline_z = make_splrep(range(len(z)), z, s=0)

            # Calculate splined coordinates
            step_values = np.linspace(
                0,
                len(carbon_alpha_coords) - 1,
                round(len(carbon_alpha_coords) * smoothness),
            )
            helix_x = splev(step_values, spline_x)
            helix_y = splev(step_values, spline_y)
            helix_z = splev(step_values, spline_z)
            helix_coords = list(zip(helix_x, helix_y, helix_z))
            coords.extend(helix_coords)

            # Colors
            hash_value = int(md5(chain.get_id().encode("utf-8")).hexdigest(), 16)
            r = (hash_value >> 16) & 0xFF
            g = (hash_value >> 8) & 0xFF
            b = hash_value & 0xFF
            chain_color = color.rgb(r / 255, g / 255, b / 255)
            colors.extend([chain_color for _ in helix_coords])

            # Segments (triangles)
            segments.extend(
                [
                    (i, i + 1)
                    for i in range(segment_start, segment_start + len(helix_coords) - 1)
                ]
            )

        return Mesh(
            mode="line",
            vertices=coords,
            colors=colors,
            triangles=segments,
            thickness=thickness,
        )
