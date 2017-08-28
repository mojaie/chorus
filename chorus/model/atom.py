#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import os
import yaml
from collections import Counter

""" Import periodic table """
with open(os.path.join(
        os.path.dirname(__file__), "periodictable.yaml")) as file:
    P_TAB = yaml.load(file.read())


def atom_number(symbol):
    return P_TAB[symbol]['number']


ALIASES = {
    "name": "name",
    "n": "number",
    "sym": "symbol",
    "chg": "charge",
    "hs": "H_count",
    "vis": "visible",
    "crds": "coords",
    "c": "color",
    "pi": "pi",
    "arom": "aromatic",
    "hdo": "H_donor",
    "hac": "H_acceptor",
    "cbnl": "carbonyl_C",
    "lp": "lone_pair",
    "wc": "wctype",
    "patty": "patty",
    "stereo": "stereo"
}


class Atom(object):
    def __init__(self, symbol, attr_dict=None):
        self.name = P_TAB[symbol]['name']
        self.number = P_TAB[symbol]['number']
        self.symbol = symbol
        self.charge = 0
        self.H_count = 0
        self.multi = 1
        self.mass = None

        # Appearance
        self.visible = 1
        if symbol == 'C':
            self.visible = 0
        self.coords = None
        self.color = tuple(P_TAB[symbol].get('color', [0, 192, 192]))

        # chem.descriptor
        self.pi = 0
        self.aromatic = 0
        self.H_donor = 0
        self.H_acceptor = 0
        self.carbonyl_C = 0
        self.lone_pair = 0  # not used?
        if self.symbol in ('N', 'O'):
            self.lone_pair = 1
        if self.symbol in ('N', 'O', 'F'):
            self.H_acceptor = 1
        self.wctype = None  # chem.wclogp
        self.patty = 7  # chem.patty
        # 1:cation 2:anion 3:donor 4:acceptor
        # 5:polar 6:hydrophobe 7:others

        # for SMILES
        self.stereo = 0
        # self.excited = False
        # self.stereo_flag = None

        # loads attribute dict
        if attr_dict is not None:
            for k, v in attr_dict.items():
                self.__setattr__(k, v)

    def __setattr__(self, name, value):
        object.__setattr__(self, ALIASES.get(name, name), value)

    def __getattr__(self, name):
        return object.__getattribute__(self, ALIASES.get(name, name))

    def add_hydrogen(self, num):
        self.H_count = num
        if num > 0 and self.symbol in ("N", "O"):
            self.H_donor = 1
        else:
            self.H_donor = 0

    def formula_html(self, reversed_=False):
        if self.H_count == 1:
            text = "H"
        elif self.H_count > 1:
            text = "H<sub>{}</sub>".format(self.H_count)
        else:
            text = ""
        seq = [self.symbol, text, self.charge_sign_html()]
        if reversed_:
            seq = reversed(seq)
        return "".join(seq)

    def composition(self):
        return Counter({self.symbol: 1, 'H': self.H_count})

    def mw(self):
        m = P_TAB[self.symbol].get('std_weight')
        mh = P_TAB['H'].get('std_weight')
        return m + mh * self.H_count

    def charge_sign(self):
        if self.charge > 0:
            sign = "+"
        elif self.charge < 0:
            sign = "–"  # en dash, not hyphen-minus
        else:
            return ""
        ab = abs(self.charge)
        if ab > 1:
            return str(ab) + sign
        return sign

    def charge_sign_html(self):
        if self.charge:
            return "<sup>{}</sup>".format(self.charge_sign())
        return ""

    def dumps(self):
        res = {}
        for k in ALIASES.keys():
            res[k] = getattr(self, k)
        return res
