import _setup

from mesh import RectRectGrid

rect = RectRectGrid(5,3)

rect.mesh((10,6))

rect.view(vertices=False,centers=False)