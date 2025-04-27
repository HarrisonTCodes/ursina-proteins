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
        helices_thickness: float = 12,
        coils_thickness: float = 3,
        chains_smoothness: float = 3,
        chain_id_color_map: dict[str, Color] = dict(),
        atom_element_color_map: dict[str, Color] = dict(),
        *args,
        **kwargs,
    ):
        parser = PDBParser()
        self.structure = parser.get_structure("protein", pdb_filepath)
        self.helices = self.get_helices(pdb_filepath)
        structure_center_of_mass = self.structure.center_of_mass()

        self.atoms_entity = Entity(
            model=self.compute_atoms_mesh(atom_element_color_map),
            origin=structure_center_of_mass,
            *args,
            **kwargs,
        )

        chain_meshes = self.compute_helices_and_coils_meshes(
            chain_id_color_map, chains_smoothness, helices_thickness, coils_thickness
        )
        self.helices_entity = Entity(
            model=chain_meshes[0],
            origin=structure_center_of_mass,
            *args,
            **kwargs,
        )
        self.coils_entity = Entity(
            model=chain_meshes[1],
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

    def compute_helices_and_coils_meshes(
        self,
        id_color_map: dict[str, Color],
        smoothness: float,
        helices_thickness: float,
        coils_thickness: float,
    ) -> list[Mesh]:
        verts = {"helices": [], "coils": []}
        tris = {"helices": [], "coils": []}
        colors = {"helices": [], "coils": []}

        for chain in self.structure.get_chains():
            carbon_alpha_coords = [
                atom.coord for atom in chain.get_atoms() if atom.get_id() == "CA"
            ]
            chain_id = chain.get_id()
            chain_segments = fill_segments(
                self.helices[chain_id], len(carbon_alpha_coords), "helices", "coils"
            )

            for segment_type, segments in chain_segments.items():
                for start, end in segments:
                    atoms = carbon_alpha_coords[start : end + 1]
                    if not atoms:
                        continue
                    tris_start = len(verts[segment_type])

                    # Vertices
                    x, y, z = zip(*atoms)
                    splines = [
                        make_splrep(
                            range(len(values)), values, s=0, k=min(3, len(values) - 1)
                        )
                        for values in [x, y, z]
                    ]

                    # Calculate splined coordinates
                    step_values = np.linspace(
                        0,
                        len(atoms) - 1,
                        round(len(atoms) * smoothness),
                    )
                    smoothed_xyz = [splev(step_values, spline) for spline in splines]
                    smoothed_coords = list(zip(*smoothed_xyz))
                    verts[segment_type].extend(smoothed_coords)

                    # Colors
                    chain_id = chain.get_id()
                    chain_color = id_color_map.get(
                        chain_id,
                        Protein.CHAIN_COLORS.get(
                            chain_id, Protein.color_from_id(chain_id)
                        ),
                    )
                    colors[segment_type].extend([chain_color for _ in smoothed_coords])

                    # Triangles
                    tris[segment_type].extend(
                        [
                            (i, i + 1)
                            for i in range(
                                tris_start, tris_start + len(smoothed_coords) - 1
                            )
                        ]
                    )

        return [
            Mesh(
                mode="line",
                vertices=verts[segment_type],
                triangles=tris[segment_type],
                colors=colors[segment_type],
                thickness=thickness,
            )
            for thickness, segment_type in zip(
                (helices_thickness, coils_thickness), ("helices", "coils")
            )
        ]

    def get_helices(self, pdb_filepath: str) -> dict[str, list[tuple[int]]]:
        helices = dict()

        with open(pdb_filepath, "r") as pdb_file:
            for line in pdb_file:
                if line.startswith("HELIX"):
                    chain_id = line[19].strip()
                    start_residue = int(line[21:25].strip()) - 1
                    end_residue = int(line[33:37].strip()) - 1
                    helix_class = int(line[38:40].strip())
                    if helix_class != 1:
                        continue

                    if chain_id in helices:
                        helices[chain_id].append((start_residue, end_residue))
                    else:
                        helices[chain_id] = [(start_residue, end_residue)]

        return helices

    @staticmethod
    def color_from_id(id: str) -> Color:
        hash_value = int(md5(id.encode("utf-8")).hexdigest(), 16)
        r = (hash_value >> 16) & 0xFF
        g = (hash_value >> 8) & 0xFF
        b = hash_value & 0xFF
        return color.rgb(r / 255, g / 255, b / 255)


def fill_segments(
    segments: list[tuple[int]], size: int, in_segment_label: str, out_segment_label: str
) -> dict[str, list[tuple[int]]]:
    segments = sorted(segments)
    result = {in_segment_label: [], out_segment_label: []}
    current = 0

    for start, end in segments:
        if current < start:
            result[out_segment_label].append((current, start))
        result[in_segment_label].append((start, end))
        current = end

    if current <= size:
        result[out_segment_label].append((current, size))

    return result
