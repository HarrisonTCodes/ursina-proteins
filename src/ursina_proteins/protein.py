from os import path

from Bio.PDB import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from ursina import Color, Vec3

from ursina_proteins.atoms import AtomsEntity
from ursina_proteins.chains import ChainsEntity


class Protein:
    """
    A class to represent a protein structure and render it as entities in Ursina.

    Attributes:
        structure: The parsed protein structure.
        helices: Dictionary mapping chain IDs to lists of helix segments.
        atoms_entity: Entity representing atoms.
        chains_entity: Entity representing helices and subentity representing coils.
        entities: List of all structural entities (atoms, helices, coils).

    Class Attributes:
        PDB_CHAIN_ID_INDEX: Index of chain ID in PDB file line.
        PDB_START_RESIDUE_INDICES: Indices of start-of-helix residue in PDB file line.
        PDB_END_RESIDUE_INDICES: Indices of end-of-helix residue in PDB file line.
    """

    PDB_CHAIN_ID_INDEX = 19
    PDB_START_RESIDUE_INDICES = (21, 25)
    PDB_END_RESIDUE_INDICES = (33, 37)

    def __init__(
        self,
        protein_filepath: str,
        legacy_pdb: bool = True,
        parser_quiet: bool = True,
        chains_smoothness: float = 3,
        helices_thickness: float = 4,
        coils_thickness: float = 1,
        chain_id_color_map: dict[str, Color] = None,
        compute_atoms: bool = True,
        atom_size: float = 0.1,
        atom_vertices: list[Vec3] = None,
        atom_triangles: list[tuple[int]] = None,
        atom_normals: list[Vec3] = None,
        atom_element_color_map: dict[str, Color] = None,
        *args,
        **kwargs,
    ):
        """
        Parse a protein structure file and create visualisation entities

        Args:
            protein_filepath: Path to the protein file.
            legacy_pdb: Whether protein is in legacy PDB (or newer mmCIF) format (default: True).
            parser_quiet: Flag to enable/disable logging on parser (default: True).
            chains_smoothness: Smoothness factor for chain rendering (default: 3).
            helices_thickness: Thickness of helix meshes (default: 4).
            coils_thickness: Thickness of coil meshes (default: 1).
            chain_id_color_map: Color mapping for chains (default: None).
            compute_atoms: Flag to enable/disable atoms computation (default: True).
            atom_size: Size of individual atoms in the atoms mesh (default: 0.1).
            atom_vertices: Base vertices to use for atom geometry (default: None).
            atom_triangles: Base triangles to use for atom geometry (default: None).
            atom_normals: Base normals to use for atom geometry (default: None).
            atom_element_color_map: Color mapping for atoms (default: None).
            *args: Arguments passed to constructor for each entity.
            **kwargs: Keyword arguments passed to constructor for each entity.
        """

        if not path.isfile(protein_filepath):
            raise FileNotFoundError(f"Protein file not found: {protein_filepath}")

        # Parse structure
        parser = PDBParser() if legacy_pdb else MMCIFParser()
        parser.QUIET = parser_quiet
        self.structure = parser.get_structure("protein", protein_filepath)
        self.helices = (
            self.get_pdb_helices(protein_filepath)
            if legacy_pdb
            else self.get_cif_helices(protein_filepath)
        )
        structure_center_of_mass = self.structure.center_of_mass()

        # Create entities
        self.entities = []

        if compute_atoms:
            self.atoms_entity = AtomsEntity(
                self.structure.get_atoms(),
                atom_size,
                atom_vertices,
                atom_triangles,
                atom_normals,
                atom_element_color_map,
                origin=Vec3(*structure_center_of_mass),
                *args,
                **kwargs,
            )
            self.entities.append(self.atoms_entity)

        self.chains_entity = ChainsEntity(
            self.structure.get_chains(),
            self.helices,
            chains_smoothness,
            helices_thickness,
            coils_thickness,
            chain_id_color_map,
            origin=Vec3(*structure_center_of_mass),
            *args,
            **kwargs,
        )
        self.entities.append(self.chains_entity)
        self.entities.append(self.chains_entity.coils_entity)

    def get_pdb_helices(self, protein_filepath: str) -> dict[str, list[tuple[int]]]:
        """
        Extract helix information for a protein from a PDB file.

        This method parses the HELIX records in a PDB file to identify
        the start and end residues of helices for each chain.

        Args:
            protein_filepath: Path to the PDB file.

        Returns:
            A dictionary mapping chain IDs to lists of helices,
            where each segment is represented as a tuple of start/end indices.
        """

        helices = dict()

        with open(protein_filepath, "r") as pdb_file:
            for line in pdb_file:
                if line.startswith("HELIX"):
                    chain_id = line[Protein.PDB_CHAIN_ID_INDEX]
                    start_residue = int(
                        line[
                            Protein.PDB_START_RESIDUE_INDICES[
                                0
                            ] : Protein.PDB_START_RESIDUE_INDICES[1]
                        ].strip()
                    )
                    end_residue = int(
                        line[
                            Protein.PDB_END_RESIDUE_INDICES[
                                0
                            ] : Protein.PDB_END_RESIDUE_INDICES[1]
                        ].strip()
                    )

                    if chain_id in helices:
                        helices[chain_id].append((start_residue, end_residue))
                    else:
                        helices[chain_id] = [(start_residue, end_residue)]

        return helices

    def get_cif_helices(
        self, protein_filepath: str
    ) -> dict[str, list[tuple[int, int]]]:
        helices = {}

        # Load CIF file into memory
        with open(protein_filepath, "r") as file:
            lines = [line.strip() for line in file if line.strip()]

        loop_start = None
        for i, line in enumerate(lines):
            # Heuristically check loop is data block
            if line.startswith("loop_") and any(
                "_struct_conf." in line for line in lines[i + 1 : i + 10]
            ):
                loop_start = i
                break
        if loop_start is None:
            return helices

        identifiers = []
        data_start = loop_start + 1
        for j, line in enumerate(lines[data_start:], start=data_start):
            # Get index where identifiers end and data begins
            if not line.startswith("_struct_conf."):
                data_start = j
                break
            # Get identifiers in order they appear in line
            identifiers.append(line.split(".")[1])

        # Create identifier-index map for more readable/semantic access and confirm key identifiers are present
        identifier_index_map = {
            identifier: k for k, identifier in enumerate(identifiers)
        }
        if (
            not {
                "conf_type_id",
                "beg_auth_asym_id",
                "beg_auth_seq_id",
                "end_auth_seq_id",
            }
            <= identifier_index_map.keys()
        ):
            return helices

        for line in lines[data_start:]:
            # End processing when next loop reached
            if line.startswith("loop_"):
                break

            parts = line.split()
            # Skip irrelevant lines
            if len(parts) < len(identifiers):
                continue
            if "HELX" not in parts[identifier_index_map["conf_type_id"]]:
                continue

            chain_id = parts[identifier_index_map["beg_auth_asym_id"]]
            start = int(parts[identifier_index_map["beg_auth_seq_id"]])
            end = int(parts[identifier_index_map["end_auth_seq_id"]])
            helices.setdefault(chain_id, []).append((start, end))

        return helices
