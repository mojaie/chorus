#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

ALIASES = {
    "o": "order",
    "vis": "visible",
    "arom": "aromatic",
    "rot": "rotatable",
    "ilf": "is_lower_first",
    "type": "type",
    "sgeo": "smiles_cis_trans"
}


class Bond(object):
    def __init__(self, attr_dict=None):
        self.visible = 1
        # bond order(1-3)  9: aromatic(for SMILES parser)
        self.order = 1
        self.aromatic = 0
        self.rotatable = 0

        """If true, lower index atom is the first atom and the other is
        second atom when the bond has directional feature.
        (If false, higher index atom is the first)
        """
        self.is_lower_first = 1

        """Bond type:
        Single bond:
            0: first - second
            1: first ◀ second (Up-arrow)
            2: first ◁ second (Down-arrow)
            3: first ~ second (Chiral)
        Double bond:
            0: second ニ first (clockwise, default)
            1: first ニ second (counter-clockwise)
            2: first ＝ second (equal length, for only terminal bond by default)
            3: first × second (Cis-Trans Unknown)
        """
        self.type = 0

        self.smiles_cis_trans = 0

        # loads attribute dict
        if attr_dict is not None:
            for k, v in attr_dict.items():
                setattr(self, k, v)

    def __setattr__(self, name, value):
        object.__setattr__(self, ALIASES.get(name, name), value)

    def __getattr__(self, name):
        return object.__getattribute__(self, ALIASES.get(name, name))

    def dumps(self):
        res = {}
        for k in ALIASES.keys():
            res[k] = getattr(self, k)
        return res
