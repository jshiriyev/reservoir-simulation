from ._items import Slot

class Survey():
    """It is a well survey (direction or trajectory)."""
    def __init__(self,**kwargs):

        self.slot = Slot(**kwargs)

        self.tops = [] # It must be a BinarySearchTree

    def add_top(self,name,depth,**kwargs):

        zone = Zone(name=name,depth=depth,**kwargs)

        self.tops.insert(zone)