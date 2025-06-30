from os import path

from ursina import Button, DirectionalLight, EditorCamera, Ursina, Vec2, window
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

    toggle_atoms_button = Button(
        text="Toggle atoms visibility",
        on_click=lambda: protein.atoms_entity.disable()
        if protein.atoms_entity.enabled
        else protein.atoms_entity.enable(),
        scale=(0.3, 0.05),
        position=window.top_left + Vec2(0.16, -0.04),
    )

    EditorCamera()
    app.run()


if __name__ == "__main__":
    main()
