import numpy as np
from collections import defaultdict
from decimal import *
getcontext().prec = 32

T_MAP = {
    'A' : np.asarray([0, 0]),
    'C' : np.asarray([0, 0.5]),
    'G' : np.asarray([0.5, 0]),
    'T' : np.asarray([0.5, 0.5]),
}

class Node(object):

    def __init__(self, char, counter = {}):
        self.char = char
        self.counter = defaultdict(Decimal, counter)
        self.entropy = {}
        self.parent = None
        self.children = {}

    def __getitem__(self, char):
        # defaultdict insert if not present
        if char not in self.children:
            self[char] = Node(char)
        return self.children[char]

    def __setitem__(self, char, node):
        self.children[char] = node
        node.parent = self

class KmerTrie(object):

    def __init__(self, k, alphabet=['A', 'T', 'C', 'G'], unk_char='?'):
        self.root = Node("")
        self.k = k
        self.alphabet = alphabet
        self.unk_char = unk_char

    def increment(self, kmer, name, value=1):
        value = Decimal(value)
        self.root.counter[name] += value
        nlist = [self.root]
        for char in kmer[::-1]:
            new_nlist = []
            for node in nlist:
                pchars = []

                # split value between all chars equally
                if char == self.unk_char:
                    value /= len(self.alphabet)
                    for b in self.alphabet:
                        pchars.append(b)
                else:
                    pchars.append(char)

                for pchar in pchars:
                    cnode = node[pchar]
                    cnode.counter[name] = cnode.counter[name] + value

                    # delete empty nodes
                    if all([v == 0 for v in cnode.counter.values()]):
                        node.children.pop(pchar)
                    new_nlist.append(cnode)
            nlist = new_nlist

    def __getitem__(self, kmer, name):
        node = self.root
        for char in kmer[::-1]:
            if char not in node:
                return None
            node = node[char]
        return node.counter.get(name, None)

    def __copy__(self):
        cls = self.__class__
        res = cls.__new__(cls)
        res.k = self.k
        res.alphabet = self.alphabet
        res.unk_char = self.unk_char
        res.root = Node("")

        # DFS copy nodes
        nlist = [(self.root, res.root)]
        while len(nlist) != 0:
            node, copy_node = nlist.pop()
            for b, bnode in node.children.items():
                copy_node[b] = Node(b, bnode.counter)
                nlist.append((bnode, copy_node[b]))

        return res

    def insert(self, seq, name, value=1):
        # insert all kmers into trie
        for i in range(self.k, len(seq) + 1):
            kmer = seq[max(0, i - self.k):i]
            self.increment(kmer, name, value)

    def delete(self, seq, name, value=1):
        self.insert(seq, name, -value)

    def relative_entropy(self, pname, qname, phi, shannon=False, all_ks=False):
        phi = Decimal(phi)
        a_size = Decimal(len(self.alphabet))
        weight = sum([phi**j for j in range(self.k + 1)])

        total_entropy = 0
        self.root.entropy[pname] = 0
        self.root.entropy[qname] = 0

        # (node, depth)
        nlist = [(n, 1) for n in self.root.children.values()]

        while len(nlist) != 0:
            node, d = nlist.pop()
            if pname in node.counter:
                node.entropy[pname] = node.parent.entropy[pname] + (a_size * phi)**d * node.counter[pname]
            else:
                node.entropy[pname] = node.parent.entropy[pname]
            if qname in node.counter:
                node.entropy[qname] = node.parent.entropy[qname] + (a_size * phi)**d * node.counter[qname]
            else:
                node.entropy[qname] = node.parent.entropy[qname]

            pdensity = (1 + (node.entropy[pname] / self.root.counter[pname])) / weight
            qdensity = (1 + (node.entropy[qname] / self.root.counter[qname])) / weight
            if d == self.k:
                if shannon:
                    total_entropy += (pdensity * Decimal(np.log(float(pdensity / qdensity)))) / a_size ** self.k
                else:
                    total_entropy += (pdensity ** 2 / qdensity) / a_size ** self.k
                continue

            for b in self.alphabet:

                # add entropy of all remaining suffixes for remaining k-d levels
                if b not in node.children or \
                    (pname not in node.children[b].counter and \
                     qname not in node.children[b].counter):
                    if shannon:
                        total_entropy += (pdensity * Decimal(np.log(float(pdensity / qdensity)))) / a_size**(d + 1)
                    else:
                        total_entropy += (pdensity ** 2 / qdensity) / a_size ** (d + 1)
                else:
                    nlist.append((node[b], d + 1))
        if shannon:
            return float(total_entropy)
        else:
            return np.log(float(total_entropy))

    def get_entropy(self, pname, phis, all_ks=False):
        res = {}
        for phi in phis:
            phi = Decimal(phi)
            a_size = Decimal(len(self.alphabet))
            weight = sum([phi**j for j in range(self.k + 1)])

            if all_ks:
                total_entropy = [0 for _ in range(self.k)]
            else:
                total_entropy = 0
            self.root.entropy = 0
            nlist = [(n, 1) for n in self.root.children.values()]

            # DFS through trie
            while len(nlist) != 0:
                node, d = nlist.pop()
                node.entropy = node.parent.entropy + (a_size * phi)**d * node.counter[pname]

                if all_ks:
                    total_entropy[d - 1] += \
                        ((1 + node.entropy / self.root.counter[pname]) / \
                            (sum([phi**j for j in range(d + 1)]))
                        )**2 / a_size**d
                    if d == self.k:
                        continue
                    for b in self.alphabet:
                        if b not in node.children or b not in node.children[b].counter:
                            for rd in range(d, self.k):
                                total_entropy[rd] += \
                                    ((1 + node.entropy / self.root.counter[pname]) /
                                        (sum([phi**j for j in range(rd + 2)]))
                                    )**2 / a_size**(rd + 1)
                        else:
                            nlist.append((node[b], d + 1))
                else:
                    # exit at k-depth leaf node
                    if d == self.k:
                        total_entropy += ((1 + node.entropy / self.root.counter[pname]) / (weight))**2 / a_size**self.k
                        continue
                    for b in self.alphabet:

                        # add entropy of all remaining suffixes for remaining k-d levels
                        if b not in node.children or b not in node.children[b].counter:
                            total_entropy += ((1 + node.entropy / self.root.counter[pname]) / weight)**2 / a_size**(d + 1)
                        else:
                            nlist.append((node[b], d + 1))
            if all_ks:
                res[float(phi)] = [-np.log(float(e)) for e in total_entropy]
            else:
                res[float(phi)] = -np.log(float(total_entropy))
        return res

    def get_density(self, phi):
        xs = []
        ys = []
        heights = []

        step_size = 1 / np.log2(len(self.alphabet))**self.k
        phi = Decimal(phi)
        a_size = Decimal(len(self.alphabet))
        weight = sum([phi**j for j in range(self.k + 1)])

        self.root.entropy = 0
        nlist = [(n, 1, T_MAP[b]) for b, n in self.root.children.items()]

        while len(nlist) != 0:
            node, d, bl = nlist.pop()
            block_size = 1 / np.log2(len(self.alphabet))**d
            node.entropy = node.parent.entropy + (a_size * phi)**d * node.counter

            # exit at k-depth leaf node
            if d == self.k:
                xs.append(bl[0])
                ys.append(bl[1])
                heights.append(float((1 + node.entropy / self.root.counter) / (weight * a_size**self.k)))
                continue
            for b in self.alphabet:

                # add entropy of all remaining suffixes for remaining k-d levels
                if b not in node.children:
                    new_bl = bl + T_MAP[b] * block_size
                    _x = np.arange(new_bl[0], new_bl[0] + 0.5 * block_size, step_size)
                    _y = np.arange(new_bl[1], new_bl[1] + 0.5 * block_size, step_size)
                    _xx, _yy = np.meshgrid(_x, _y)
                    _xx = _xx.ravel()
                    _yy = _yy.ravel()
                    xs.extend(list(_xx))
                    ys.extend(list(_yy))
                    heights.extend([float((1 + node.entropy / self.root.counter) / (weight * a_size**self.k))] * len(_xx))
                else:
                    nlist.append((node[b], d + 1, bl + T_MAP[b] * block_size))

        return xs, ys, heights
