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
    """Bond object

    * Bond type
        * Single bond
            * 0: first - second
            * 1: first ◀ second (Up-arrow)
            * 2: first ◁ second (Down-arrow)
            * 3: first ~ second (Chiral)

        *Double bond
            * 0: second ニ first (clockwise, default)
            * 1: first ニ second (counter-clockwise)
            * 2: first ＝ second (equal length, for terminal bond by default)
            * 3: first × second (Cis-Trans Unknown)

    Parameters:
        attr_dict (dict): load ``Bond`` object from dict type notation.

    Attributes:
        order (int): bond order. ``smilesparser`` assigns order=9 to
            aromatic bonds.
        is_lower_first (bool): bond direction. If true, lower index atom is
            the first atom and the other is the second atom. ``Bond.type``
            attribute uses this feature for stereochemistry determination.
        type: bond type (see above)
        rotatable (bool): rotatable or not
        aromatic (bool): aromatic or not
        smiles_cis_trans: SMILES cis/trans flag
        visible (bool): whether it should be drawn or not.


    """
    def __init__(self, attr_dict=None):
        self.order = 1
        self.is_lower_first = 1
        self.type = 0
        self.rotatable = 0
        self.aromatic = 0
        self.smiles_cis_trans = 0
        self.visible = 1

        # loads attribute dict
        if attr_dict is not None:
            for k, v in attr_dict.items():
                setattr(self, ALIASES[k], v)

    def __getattr__(self, name):
        if name not in ALIASES:
            raise AttributeError()
        return object.__getattribute__(self, ALIASES[name])

    def dumps(self):
        res = {}
        for k, v in ALIASES.items():
            res[k] = getattr(self, v)
        return res
