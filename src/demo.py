from ursina import EditorCamera, Ursina

from ursina_proteins.protein import Protein


def main():
    pdb_filepath = input("Enter PDB file path: ")
    app = Ursina(borderless=False)

    Protein(pdb_filepath)

    EditorCamera()
    app.run()


if __name__ == "__main__":
    main()
