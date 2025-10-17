from math import sqrt
from typing import Iterator

from Bio.PDB.Atom import Atom
from ursina import Color, Entity, Mesh, Vec3, color

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

ELEMENT_COLORS = {
    "H": color.rgb(0.8, 0.8, 0.8),
    "C": color.rgb(0.2, 0.2, 0.2),
    "N": color.rgb(0, 0, 0.8),
    "O": color.rgb(0.8, 0, 0),
    "S": color.rgb(0.8, 0.8, 0),
    "P": color.rgb(1, 0.65, 0),
    "Cl": color.rgb(0, 0.8, 0),
    "Fe": color.rgb(0.7, 0.45, 0.2),
}


class AtomsEntity(Entity):
    """
    A class to represent atoms as a Ursina Entity.
    """

    def __init__(
        self,
        atoms: Iterator[Atom],
        atom_size: float,
        atom_vertices: list[Vec3] = None,
        atom_triangles: list[tuple[int]] = None,
        atom_normals: list[Vec3] = None,
        element_color_map: dict[str, Color] = None,
        *args,
        **kwargs,
    ):
        """
        Create an entity to visualise atoms.

        Args:
            atoms: Atoms to visualise.
            atom_size: Size of individual atoms.
            atom_vertices: Base vertices to use for atom geometry (default: None).
            atom_triangles: Base triangles to use for atom geometry (default: None).
            atom_normals: Base normals to use for atom geometry (default: None).
            element_color_map: Color mapping for atom elements (default: None).
            *args: Arguments passed to constructor for the entity.
            **kwargs: Keyword arguments passed to constructor for the entity.
        """
        if atom_size <= 0:
            raise ValueError("Atom size must be positive")

        atom_vertices = atom_vertices or ICOSAHEDRON_VERTS
        atom_triangles = atom_triangles or ICOSAHEDRON_FACES
        atom_normals = atom_normals or ICOSAHEDRON_NORMALS
        element_color_map = element_color_map or dict()

        verts = []
        faces = []
        colors = []
        norms = []

        for index, atom in enumerate(atoms):
            # Vertices
            verts.extend(
                [(vert * atom_size) + atom.get_coord() for vert in atom_vertices]
            )

            # Faces (triangles)
            faces.extend(
                [
                    tuple(i + len(atom_vertices) * index for i in face)
                    for face in atom_triangles
                ]
            )

            # Colors
            colors.extend(
                [
                    element_color_map.get(
                        atom.element,
                        ELEMENT_COLORS.get(atom.element, color.rgb(1, 0.7, 0.8)),
                    )
                    for _ in atom_vertices
                ]
            )

            # Normals
            norms.extend(atom_normals)

        atoms_mesh = Mesh(vertices=verts, triangles=faces, colors=colors, normals=norms)
        super().__init__(model=atoms_mesh, *args, **kwargs)
