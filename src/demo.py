from ursina import EditorCamera, Ursina

from ursina_proteins.protein import Protein


def main():
    app = Ursina()

    Protein()

    EditorCamera()
    app.run()


if __name__ == "__main__":
    main()
