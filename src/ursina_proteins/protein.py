from hashlib import md5
from math import sqrt

import numpy as np
from Bio.PDB import PDBParser
from scipy.interpolate import make_splrep, splev
from ursina import Color, Entity, Mesh, Vec3, color

# Geometry constants
PHI = (1 + sqrt(5)) / 2

ICOSAHEDRON_VERTS = [
    Vec3(-1, PHI, 0),
    Vec3(1, PHI, 0),
    Vec3(-1, -PHI, 0),
    Vec3(1, -PHI, 0),
    Vec3(0, -1, PHI),
    Vec3(0, 1, PHI),
    Vec3(0, -1, -PHI),
    Vec3(0, 1, -PHI),
    Vec3(PHI, 0, -1),
    Vec3(PHI, 0, 1),
    Vec3(-PHI, 0, -1),
    Vec3(-PHI, 0, 1),
]

ICOSAHEDRON_FACES = [
    (0, 11, 5),
    (0, 5, 1),
    (0, 1, 7),
    (0, 7, 10),
    (0, 10, 11),
    (1, 5, 9),
    (5, 11, 4),
    (11, 10, 2),
    (10, 7, 6),
    (7, 1, 8),
    (3, 9, 4),
    (3, 4, 2),
    (3, 2, 6),
    (3, 6, 8),
    (3, 8, 9),
    (4, 9, 5),
    (2, 4, 11),
    (6, 2, 10),
    (8, 6, 7),
    (9, 8, 1),
]

ICOSAHEDRON_NORMALS = [v.normalized() for v in ICOSAHEDRON_VERTS]


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

    CHAIN_COLORS = {
        "A": color.rgb(1, 0, 0),
        "B": color.rgb(0, 1, 0),
        "C": color.rgb(0, 0, 1),
        "D": color.rgb(1, 1, 0),
        "E": color.rgb(1, 0.5, 0.8),
        "F": color.rgb(0.2, 0.7, 1),
        "G": color.rgb(1, 0.6, 0),
        "H": color.rgb(1, 0, 1),
    }

    def __init__(
        self,
        pdb_filepath: str,
        chains_thickness: float = 4,
        chains_smoothness: float = 3,
        chain_id_color_map: dict[str, Color] = dict(),
        atom_element_color_map: dict[str, Color] = dict(),
        *args,
        **kwargs,
    ):
        parser = PDBParser()
        self.structure = parser.get_structure("protein", pdb_filepath)
        structure_center_of_mass = self.structure.center_of_mass()

        self.atoms_entity = Entity(
            model=self.compute_atoms_mesh(atom_element_color_map),
            origin=structure_center_of_mass,
            *args,
            **kwargs,
        )

        self.chains_entity = Entity(
            model=self.compute_chains_mesh(
                chain_id_color_map, chains_thickness, chains_smoothness
            ),
            origin=structure_center_of_mass,
            *args,
            **kwargs,
        )

    def compute_atoms_mesh(self, element_color_map: dict[str, Color]) -> Mesh:
        verts = []
        faces = []
        colors = []
        norms = []

        for index, atom in enumerate(self.structure.get_atoms()):
            # Vertices
            verts.extend(
                [(vert * 0.1) + atom.get_coord() for vert in ICOSAHEDRON_VERTS]
            )

            # Faces (triangles)
            faces.extend(
                [
                    tuple(i + len(ICOSAHEDRON_VERTS) * index for i in face)
                    for face in ICOSAHEDRON_FACES
                ]
            )

            # Colors
            colors.extend(
                [
                    element_color_map.get(
                        atom.element,
                        Protein.ELEMENT_COLORS.get(
                            atom.element, color.rgb(1, 0.7, 0.8)
                        ),
                    )
                    for _ in ICOSAHEDRON_VERTS
                ]
            )

            # Normals
            norms.extend(ICOSAHEDRON_NORMALS)

        return Mesh(vertices=verts, triangles=faces, colors=colors, normals=norms)

    def compute_chains_mesh(
        self, id_color_map: dict[str, Color], thickness: float, smoothness: float
    ) -> Mesh:
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
            chain_id = chain.get_id()
            chain_color = id_color_map.get(
                chain_id,
                Protein.CHAIN_COLORS.get(chain_id, Protein.color_from_id(chain_id)),
            )
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

    @staticmethod
    def color_from_id(id: str) -> Color:
        hash_value = int(md5(id.encode("utf-8")).hexdigest(), 16)
        r = (hash_value >> 16) & 0xFF
        g = (hash_value >> 8) & 0xFF
        b = hash_value & 0xFF
        return color.rgb(r / 255, g / 255, b / 255)
