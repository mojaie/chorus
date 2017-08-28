#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import networkx as nx

from chorus.model.atom import Atom
from chorus.model.bond import Bond


class Compound(object):
    """molecule in graph expression.

    Attributes:
        graph (networkx.Graph): molecule graph
        data (dict): optional attributes(ex. general name, safety info)
        descriptors (set): Available descriptors
    """
    def __init__(self, data=None):
        self.graph = nx.Graph()
        self.data = {}
        self.descriptors = set()

        """ Topology (topology) """
        self.rings = None
        self.scaffolds = None
        self.isolated = None

        """ Canvas size (draw.helper) """
        self.size2d = None

        # loads JSON formatted dict
        if data is not None:
            if "data" in data:
                self.data = data["data"]
            else:  # for backword compatibility (will be removed)
                self.data = data["options"]
            self.descriptors = set(data["descriptors"])
            self.rings = data["rings"]
            self.scaffolds = data["scaffolds"]
            self.isolated = data["isolated"]
            self.size2d = data["size2d"]
            for n, atom in data["atoms"].items():
                self.add_atom(int(n), Atom(atom["sym"], atom))
            for u, conn in data["connections"].items():
                for v, bond in conn.items():
                    self.add_bond(int(u), int(v), Bond(bond))

    def __str__(self):
        """Simple molecule description"""
        adscs = []
        for i, a in self.atoms_iter():
            adsc = '{}({}{} {} {})'.format(
                i, a.symbol, a.charge_sign(), int(a.pi), a.stereo)
            adscs.append(adsc)
        atext = ', '.join(adscs)
        odsign = {1: '-', 2: '=', 3: '#'}
        bdscs = []
        for u, v, b in self.bonds_iter():
            bdsc = '{}{}{}'.format(u, odsign[b.order], v)
            bdscs.append(bdsc)
        btext = ', '.join(bdscs)
        return '\n'.join(['Compound:', atext, btext])

    def __len__(self):
        """Alias of atom_count"""
        return self.atom_count()

    def require(self, desc):
        if desc not in self.descriptors:
            raise TypeError("Descriptor '{}' is required.".format(desc))

    def atom(self, key):
        """Get an atom."""
        return self.graph.node[key]['atom']

    def add_atom(self, key, atom):
        """Set an atom. Existing atom will be overwritten."""
        self.graph.add_node(key, attr_dict={'atom': atom})

    def remove_atom(self, key):
        """Remove an atom and adjacent bonds."""
        self.graph.remove_node(key)

    def atoms_iter(self):
        """Iterate over atoms."""
        for n, attr in self.graph.nodes_iter(data=True):
            yield n, attr['atom']

    def atom_count(self):
        """Get number of atoms."""
        return self.graph.number_of_nodes()

    def bond(self, key1, key2):
        """Get a bond."""
        return self.graph.edge[key1][key2]['bond']

    def add_bond(self, key1, key2, bond):
        """Set a bond. Existing bond will be overwritten."""
        self.graph.add_edge(key1, key2, attr_dict={'bond': bond})

    def remove_bond(self, key1, key2):
        """Remove a bond."""
        self.graph.remove_edge(key1, key2)

    def bonds_iter(self):
        """Iterate over bonds."""
        for u, v, attr in self.graph.edges_iter(data=True):
            yield u, v, attr['bond']

    def bond_count(self):
        """Return number of bonds."""
        return self.graph.number_of_edges()

    def key_set(self):
        """Get a set of atom keys"""
        return set(self.graph.adj.keys())

    def neighbors(self, key):
        """Return dict of neighbor atom index and connecting bond."""
        return {n: attr['bond'] for n, attr in self.graph[key].items()}

    def neighbor_count(self, key):
        """Return number of neighbors."""
        return len(self.graph[key])

    def neighbors_iter(self):
        """Iterate over atoms and return its neighbors."""
        for n in self.graph:
            yield n, self.neighbors(n)

    def clear(self):
        """Empty the instance """
        # self.graph = nx.Graph()
        self.graph.clear()
        self.data.clear()
        self.descriptors.clear()
        self.size2d = None
        self.rings = None
        self.scaffolds = None
        self.isolated = None

    def add_molecule(self, mol, bond=None, base=None, target=None):
        """connect atom group (for SMILES parser)

        May requires recalculation of 2D coordinate for drawing

        Args:
            mol: graphmol.Compound()
                the original object will be copied.
            bond: Bond object to be connected.
                the original will not be copied so be careful.
            base: index of atom in self to connect
            target: index of atom in group to be connected
        Raises:
            TypeError
        """
        ai = self.available_idx()
        mapping = {n: n + ai - 1 for n, _ in mol.atoms_iter()}
        relabeled = nx.relabel_nodes(mol.graph, mapping)  # copy=True
        self.graph.add_nodes_from(relabeled.nodes(data=True))
        self.graph.add_edges_from(relabeled.edges(data=True))
        if bond:
            self.add_bond(base, mapping[target], bond)

    def available_idx(self):
        return max(self.graph.nodes() + [0, ]) + 1

    def expand(self, idx):
        """ expand shorthand symbol to graphmol """
        # TODO: not yet implemented
        pass

    def collapse(self, idxs, newidx, symbol):
        pass

    def jsonized(self):
        data = {
            "descriptors": list(self.descriptors),
            "isolated": self.isolated,
            "data": self.data,
            'rings': self.rings,
            'scaffolds': self.scaffolds,
            'size2d': self.size2d
        }
        data["atoms"] = {}
        data["connections"] = {}
        for n, atom in self.atoms_iter():
            data["atoms"][n] = atom.dumps()
        for u, v, bond in self.bonds_iter():
            if u not in data["connections"]:
                data["connections"][u] = {}
            data["connections"][u][v] = bond.dumps()
        return data
