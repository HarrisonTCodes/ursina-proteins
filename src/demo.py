from os import path

from ursina import DirectionalLight, EditorCamera, Ursina
from ursina.shaders import lit_with_shadows_shader

from ursina_proteins.protein import Protein


def main():
    app = Ursina(borderless=False)
    base_dir = path.dirname(__file__)

    protein = Protein(
        path.join(base_dir, "..", "assets", "insulin.pdb"),
        shader=lit_with_shadows_shader,
    )
    light = DirectionalLight(position=(1000, 1000, 1000))
    light.look_at(protein.atoms_entity)

    EditorCamera()
    app.run()


if __name__ == "__main__":
    main()
