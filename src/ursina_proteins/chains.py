from hashlib import md5
from typing import Iterator

import numpy as np
from Bio.PDB.Chain import Chain
from scipy.interpolate import make_splrep, splev
from ursina import Color, Entity, Mesh, color

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


class ChainsEntity(Entity):
    """
    A class to represent chains as a Ursina Entity.

    Attributes:
        coils_entity: Subentity representing coils.
    """

    def __init__(
        self,
        chains: Iterator[Chain],
        helices: dict[str, list[tuple[int]]],
        smoothness: float,
        helices_thickness: float,
        coils_thickness: float,
        id_color_map: dict[str, Color],
        *args,
        **kwargs,
    ):
        """
        Create an entity to visualise helices and a subentity to visualise coils of chains

        Args:
            chains: Chains to visualise.
            helices: Map of chain IDs to chain helix atom indices.
            smoothness: Factor controlling the smoothness of the chains.
            helices_thickness: Thickness of helix meshes.
            coils_thickness: Thickness of coil meshes.
            id_color_map: Color mapping for chain IDs.
            *args: Arguments passed to constructor for the entities.
            **kwargs: Keyword arguments passed to constructor for the entities.
        """

        if helices_thickness <= 0 or coils_thickness <= 0:
            raise ValueError("Thickness values must be positive")

        if smoothness < 1:
            raise ValueError("Smoothness value must be at least 1")

        verts = {"helices": [], "coils": []}
        tris = {"helices": [], "coils": []}
        colors = {"helices": [], "coils": []}

        for chain in chains:
            # Map of atom number to atom coordinate
            carbon_alpha_coords = {
                atom.get_parent().get_id()[1]: atom.coord
                for atom in chain.get_atoms()
                if atom.get_id() == "CA"
            }

            # Chain info
            chain_id = chain.get_id()
            chain_helices = helices.get(chain_id) or []
            chain_segments = self.parse_segments(
                chain_helices, len(carbon_alpha_coords)
            )

            # Render each segment (helices and coils)
            for segment_type, segments in chain_segments.items():
                for start, end in segments:
                    # Get coordinates of the segment's carbon alpha atoms
                    coords = [
                        coord
                        for i in range(start, end + 1)
                        if (coord := carbon_alpha_coords.get(i)) is not None
                    ]

                    tris_start = len(verts[segment_type])

                    # Vertices
                    x, y, z = zip(*coords)
                    splines = [
                        make_splrep(
                            range(len(values)), values, s=0, k=min(3, len(values) - 1)
                        )
                        for values in [x, y, z]
                    ]

                    # Calculate splined coordinates
                    step_values = np.linspace(
                        0,
                        len(coords) - 1,
                        round(len(coords) * smoothness),
                    )
                    smoothed_xyz = [splev(step_values, spline) for spline in splines]
                    smoothed_coords = list(zip(*smoothed_xyz))
                    verts[segment_type].extend(smoothed_coords)

                    # Colors
                    chain_color = id_color_map.get(
                        chain_id,
                        CHAIN_COLORS.get(chain_id, self.color_from_id(chain_id)),
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

        helices_mesh = Mesh(
            mode="line",
            vertices=verts["helices"],
            triangles=tris["helices"],
            colors=colors["helices"],
            thickness=helices_thickness,
        )
        coils_mesh = Mesh(
            mode="line",
            vertices=verts["coils"],
            triangles=tris["coils"],
            colors=colors["coils"],
            thickness=coils_thickness,
        )

        super().__init__(model=helices_mesh, *args, **kwargs)
        self.coils_entity = Entity(model=coils_mesh, *args, **kwargs)

    @staticmethod
    def color_from_id(id: str) -> Color:
        """
        Generate a deterministic color based on a string identifier.

        This method creates a consistent color for a given ID string by hashing
        the string and extracting RGB values from the hash.

        Args:
            id: String identifier to generate a color for.

        Returns:
            A Color object with RGB values derived from the hash of the input ID.
        """

        hash_value = int(md5(id.encode("utf-8")).hexdigest(), 16)
        r = (hash_value >> 16) & 0xFF
        g = (hash_value >> 8) & 0xFF
        b = hash_value & 0xFF
        return color.rgb(r / 255, g / 255, b / 255)

    @staticmethod
    def parse_segments(
        segments: list[tuple[int]],
        size: int,
    ) -> dict[str, list[tuple[int]]]:
        """
        Parse a list of segments and fill in the gaps between them.

        This utility function takes a list of segment indices and generates
        a complete segmentation by filling in the gaps between them.

        Args:
            segments: List of segments, each represented as a tuple of (start, end) indices.
            size: The total size to cover.

        Returns:
            A dictionary with two keys ("helices" and "coils"),
            each mapping to a list of (start, end) tuples representing the segments.
        """

        segments = sorted(segments)
        result = {"helices": [], "coils": []}
        current = 0

        for start, end in segments:
            if current < start:
                result["coils"].append((current, start))
            result["helices"].append((start, end))
            current = end

        if current <= size:
            result["coils"].append((current, size))

        return result
