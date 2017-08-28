#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

from chorus.model.Abstractatom import AbstractAtom


def collapse(self):
    """ recognize amino acid moiety and convert into AminoAcid object """
    # TODO: not implemented
    pass


class AminoAcid(AbstractAtom):
    def __init__(self):
        self.symbol = None
        self.composition
        self.mw
        self.n_term = None
        self.c_term = None

    def formula_html(self):
        return self.symbol

    def composition(self):
        return self.composition

    def mw(self):
        return self.mw


class Glycine(AminoAcid):
    def __init__(self):
        self.name = "Glycine"
        self.symbol = "Gly"
        self.shorthand = "G"
        self.composition = None
        self.mw = None
