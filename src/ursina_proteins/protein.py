from ursina import Entity


class Protein(Entity):
    def __init__(self, *args, **kwargs):
        super().__init__(model="cube", *args, **kwargs)
