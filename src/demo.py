from ursina import DirectionalLight, EditorCamera, Ursina
from ursina.shaders import lit_with_shadows_shader

from ursina_proteins.protein import Protein


def main():
    pdb_filepath = input("Enter PDB file path: ")
    app = Ursina(borderless=False)

    protein = Protein(pdb_filepath, shader=lit_with_shadows_shader)
    light = DirectionalLight(position=(1000, 1000, 1000))
    light.look_at(protein.atoms_entity)

    EditorCamera()
    app.run()


if __name__ == "__main__":
    main()
