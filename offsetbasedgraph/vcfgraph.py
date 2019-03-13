from .vcfmap import *
from collections import defaultdict, namedtuple
import numpy as np
import logging

char_codes = {"a": 0, "c": 1, "g": 2, "t": 3}
char_rev = np.array(["a", "c", "g", "t"])

Interval = namedtuple("Interval", ["start", "end", "node_ids",  "direction"])
Position = namedtuple("Position", ["node_id", "offset"])


class VCFMap:
    def __init__(self, ids, codes):
        self._ids = ids
        self._codes = codes

    @classmethod
    def empty(cls, n):
        return cls(np.zeros(n), np.zeros(n))

    def fill_snps(self, var_ids, positions):
        args = np.argsort(positions)
        ids = np.arange(args.size)[args]
        self._ids[var_ids] = ids
        self._codes[var_ids] = 1

    def fill_insertions(self, var_ids, node_ids):
        self._ids[var_ids] = node_ids
        self._codes[var_ids] = 0


class SNPs:
    def __init__(self, node_index=[0], snps=[], seqs=[], char_seqs=True):
        self._node_index = np.asanyarray(node_index, dtype="int")
        self._snps = np.asanyarray(snps, dtype="int")
        if char_seqs:
            self._seqs = np.array([char_codes[c] for c in seqs], dtype="uint8")
        else:
            self._seqs = np.asanyarray(seqs, dtype="uint8")

    def __eq__(self, other):
        if not np.all(self._node_index == other._node_index):
            return False
        return np.all(self._snps == other._snps)

    def find_snp(self, node, offset, seq):
        seq = char_codes[seq]
        a = self._node_index[node]
        if node+1 < self._node_index.size:
            b = self._node_index[node+1]
            offsets = self._snps[a:b]
        else:
            offsets = self._snps[a:]
        valid = a + np.flatnonzero(offsets == offset)
        seq_valid = [v for v in valid if self._seqs[v] == seq]
        if not seq_valid:
            return None
        assert len(seq_valid) >= 1, (node, offset, seq, seq_valid)

        return seq_valid[0]

    def save(self, basename):
        np.save(basename+"_snps_node_index.npy", self._node_index)
        np.save(basename+"_snps_snps.npy", self._snps)
        np.save(basename+"_snps_seqs.npy", self._seqs)

    @classmethod
    def load(cls, basename):
        node_index = np.load(basename+"_snps_node_index.npy")
        snps = np.load(basename+"_snps_snps.npy")
        seqs = np.load(basename+"_snps_seqs.npy")
        return cls(node_index, snps, seqs, char_seqs=False)


class Path:
    def __init__(self, node_ids, distance_to_node):
        self._node_ids = np.asanyarray(node_ids, dtype="int")
        self._distance_to_node = np.asanyarray(distance_to_node, dtype="int")

    @classmethod
    def from_nodes_and_graph(cls, nodes, graph):
        distance_to_node = np.r_[0, np.cumsum(graph._node_lens[nodes])]
        return cls(nodes, distance_to_node)

    def get_node_intervals(self):
        return zip(self._node_ids, self._distance_to_node[:-1], self._distance_to_node[1:])

    def save(self, basename):
        np.save(basename+"_node_ids.npy", self._node_ids)
        np.save(basename+"_distance_to_node.npy", self._distance_to_node)

    def distance_to_node_id(self, node_id):
        idx = np.flatnonzero(self._node_ids == node_id)[0]
        return self._distance_to_node[idx]

    def project_interval(self, interval):
        start_idx = np.search_sorted(self._distance_to_node, interval.node_ids[0], side="left")
        end_idx = np.search_sorted(self._distance_to_node, interval.node_ids[-1], side="right")
        start = self._distance_to_node[start_idx]
        end = self._distance_to_node[end_idx]
        if interval.node_ids[0] == self._node_ids[start_idx]:
            start += interval.start
        if interval.node_ids[-1] == self._node_ids[end_idx-1]:
            end = self._distance_to_node[end_idx-1]+interval.end
        return (start, end)

    @classmethod
    def load(cls, basename):
        return cls(np.load(basename+"_node_ids.npy"),
                   np.load(basename+"_distance_to_node.npy"))


class IndexedPath(Path):
    def __init__(self, node_ids, distance_to_node, n_nodes=None, node_lookup=None):
        super().__init__(node_ids, distance_to_node)
        if node_lookup is not None:
            self._node_lookup = np.asanyarray(node_lookup, dtype="int")
            return
        if n_nodes is None:
            n_nodes = np.max(self._node_ids)+1
        self._node_lookup = -1*np.ones(n_nodes, dtype="int")
        self._node_lookup[self._node_ids] = np.arange(self._node_ids.size)

    def is_in_path(self, node_id):
        return self._node_lookup[node_id] >= 0

    def project_position(self, position):
        if self.is_in_path(position.node_id):
            return self.distance_to_node_id(position.node_id)+position.offset
        idx = np.searchsorted(self._node_ids, position.node_id)
        return self._distance_to_node[idx]

    def project_positions(self, positions):
        node_ids = positions[:, 0]
        offsets = positions[:, 1]
        is_in_path = self.is_in_path(node_ids)
        path_positions = self.distance_to_node_id(node_ids)+offsets
        idxs = np.searchsorted(self._node_ids, node_ids)
        non_path_positions = self._distance_to_node[idxs]
        return np.where(is_in_path, path_positions, non_path_positions)

    def back_project_positions(self, positions):
        idxs = np.searchsorted(self._distance_to_node, positions, side="right")-1
        offsets = np.maximum(0, positions-self._distance_to_node[idxs])
        node_ids = self._node_ids[idxs]
        return node_ids, offsets

    def distance_to_node_id(self, node_id):
        idx = self._node_lookup[node_id]
        return self._distance_to_node[idx]

    def traverse_ref_nodes(self, node_id=None):
        if node_id is None:
            node_id = self._node_ids[0]
        idx = 0
        if node_id is not None:
            idx = self._node_lookup[node_id]
        return self._node_ids[idx:]
        return (node for node in self._node_ids[idx:])
        # while node_id != -1:
        #     yield node_id
        #     node_id = self._node_ids[self._node_lookup[node_id]]

    def ref_nodes_between(self, start_node, end_node):
        start_idxs = self._node_lookup[start_node]
        end_idx = self._node_lookup[end_node]
        return self._node_ids[start_idxs:end_idx]

    def next_node(self, node_id):
        return self._node_ids[self._node_lookup[node_id]+1]

    @classmethod
    def from_indexed_interval(cls, indexed_interval):
        node_ids = indexed_interval.get_sorted_nodes_in_interval()-indexed_interval.min_node
        distances = np.hstack((
            indexed_interval._node_to_distance[node_ids],
            indexed_interval._length))
        return cls(node_ids, distances)


class Sequences:
    _letters = np.array(["n", "a", "c", "t", "g", "m"])

    def __init__(self, indices, sequence):
        self._node_indices = np.asanyarray(indices, dtype="int")
        self._sequences = np.asanyarray(sequence, dtype="uint8")

    def __getitem__(self, node_id):
        if isinstance(node_id, list) or isinstance(node_id, np.ndarray):
            return [self[n] for n in node_id]
        start, end = self._node_indices[[node_id, node_id+1]]
        return "".join(self._letters[self._sequences[start:end]])

    @classmethod
    def from_sequence_graph(cls, seq_graph):
        return  cls(seq_graph._indices[1:],
                    seq_graph._sequence_array)


class AdjList:
    def __init__(self, node_index=[0], to_nodes=[]):
        self._node_index = np.asanyarray(node_index, dtype="int")
        self._to_nodes = np.asanyarray(to_nodes, dtype="int")

    def items(self):
        return zip(self.keys(), self.values())

    def values(self):
        return (self._to_nodes[self._node_index[i]:self._node_index[i+1]]
                for i in self.keys())

    def keys(self):
        return range(self._node_index.size-1)

    def __len__(self):
        return len(self._node_index)-1

    def __getitem__(self, key):
        a = self._node_index[key]
        b = self._node_index[key+1]
        return self._to_nodes[a:b]

    def __eq__(self, other):
        if not np.all(self._node_index == other._node_index):
            logging.debug("Different node index")
            return False
        if not np.all(self._to_nodes == other._to_nodes):
            logging.debug("Different to nodes")
            return False
        logging.debug("Same adj list")
        return True

    def save(self, basename):
        np.save(basename+"_adj_list_node_index.npy", self._node_index)
        np.save(basename+"_adj_list_to_nodes.npy", self._to_nodes)

    @classmethod
    def load(cls, basename):
        node_index = np.load(basename+"_adj_list_node_index.npy")
        to_nodes = np.load(basename+"_adj_list_to_nodes.npy")
        return cls(node_index, to_nodes)

    @classmethod
    def from_dict(cls, tmp_adj_list, n_nodes):
        adj_list = defaultdict(list)
        adj_list.update(tmp_adj_list)
        node_index = np.empty(n_nodes+1, dtype="int")
        lens = [len(adj_list[from_node]) if from_node in adj_list else 0 for from_node in range(n_nodes)]
        node_index[1:] = np.cumsum(lens)
        node_index[0] = 0
        to_nodes = [to_node for from_node in range(n_nodes)
                    for to_node in sorted(adj_list[from_node])]
        return cls(node_index, to_nodes)

    def _get_reverse_edges(self):
        reverse_edges = defaultdict(list)
        for block, edges in self.items():
            for edge in edges:
                reverse_edges[edge].append(block)

        return reverse_edges

    def get_reverse(self):
        return self.from_dict(self._get_reverse_edges(), self._node_index.size-1)

    def __repr__(self):
        return "AdjList(%s, %s)" % (self._node_index, self._to_nodes)

    @classmethod
    def from_obg_adj(cls, adj_list, n_nodes):
        padding = adj_list._values.size*np.ones(n_nodes+1-adj_list._indices.size)
        indices = np.hstack((adj_list._indices, padding))
        print(adj_list._indices.size, indices.size)
        to_nodes = adj_list._values-adj_list.node_id_offset
        return cls(indices, to_nodes)


class VCFGraph:
    def __init__(self, node_lens, adj_list, snps, seqs=None):
        self._node_lens = np.asanyarray(node_lens, dtype="int")
        self._adj_list = adj_list
        self._rev_adj_list = adj_list.get_reverse()
        self._snps = snps
        self._seqs = seqs  # np.asanyarray(seqs)

    def save(self, basename):
        np.save(basename+"_node_lens.npy", self._node_lens)
        self._adj_list.save(basename)
        self._snps.save(basename)
        if self._seqs is not None:
            np.save(basename+"_seqs.npy", self._seqs)

    @classmethod
    def load(cls, basename):
        node_lens = np.load(basename+"_node_lens.npy")
        adj_list = AdjList.load(basename)
        snps = SNPs.load(basename)
        seqs = np.load(basename+"_seqs.npy")
        for i in range(seqs.size):
            seqs[i] = seqs[i].lower()
        return cls(node_lens, adj_list, snps, seqs)

    def __eq__(self, other):
        if not np.all(self._node_lens == other._node_lens):
            logging.debug("Different lens")
            return False
        if not self._adj_list == other._adj_list:
            logging.debug("Different adj list")
            return False

        if not self._snps == other._snps:
            logging.debug("Different snps")
            return False
        return np.all(self._seqs == other._seqs)

    def __repr__(self):
        if self._node_lens.size > 100:
            return "VCFGraph(%s, %s, %s)" % (self._node_lens.size, self._adj_list._to_nodes.size, sum(len(seq) for seq in self._seqs))
        return "VCFGraph(%s, %s, %s)" % (self._node_lens, self._adj_list, self._seqs)


    @classmethod
    def from_obg_graph(cls, graph, seq_graph=None):
        node_lens = graph.blocks._array[1:]
        n_nodes = node_lens.size
        adj_list = AdjList.from_obg_adj(graph.adj_list, n_nodes)
        assert adj_list._node_index.size == n_nodes+1, (adj_list._node_index.size, n_nodes)

        if seq_graph is not None:
            seqs = Sequences.from_sequence_graph(seq_graph)
            assert seqs._node_indices.size == n_nodes+1, (seqs._node_indices.size, n_nodes)

        else:
            seqs = None
        
        return cls(node_lens, adj_list, snps=SNPs(), seqs=seqs)


def classify_vcf_entry(vcf_entry):
    if len(vcf_entry.alt) > len(vcf_entry.ref):
        return INS
    elif len(vcf_entry.alt) < len(vcf_entry.ref):
        return DEL
    return SNP


def get_node_starts(deletions, insertions):
    logging.info("Getting node starts (%s del, %s ins)", len(deletions[0]), len(insertions[0]))
    break_points = np.concatenate((deletions[0],
                                   deletions[0]+deletions[1],
                                   insertions[0]))
    args = np.argsort(break_points, kind="mergesort")
    sorted_break_points = break_points[args]
    if not sorted_break_points.size:
        return np.array([0])
    diffs = np.diff(sorted_break_points)
    node_starts = np.r_[0, sorted_break_points[:-1][diffs > 0],
                        sorted_break_points[-1]]
    return node_starts


def create_adj_list(reference_node_ids, node_starts, insertions, insertion_node_ids, deletions):
    logging.info("Creating adj list")
    adj_list = defaultdict(list)

    # Add reference neighbours
    for from_node, to_node in zip(reference_node_ids[:-1], reference_node_ids[1:]):
        adj_list[from_node].append(to_node)

    node_start_idxs = np.searchsorted(node_starts, insertions[0])
    from_nodes = reference_node_ids[node_start_idxs-1]
    to_nodes = reference_node_ids[node_start_idxs]
    for from_node, to_node, node_id in zip(from_nodes, to_nodes, insertion_node_ids):
        adj_list[from_node].append(node_id)
        adj_list[node_id].append(to_node)

    node_start_idxs = np.searchsorted(node_starts, deletions[0])
    from_nodes = reference_node_ids[node_start_idxs-1]
    node_start_idxs = np.searchsorted(node_starts, deletions[0]+deletions[1])
    to_nodes = reference_node_ids[node_start_idxs]
    for from_node, to_node in zip(from_nodes, to_nodes):
        adj_list[from_node].append(to_node)
    return adj_list


def create_node_lens(node_starts, reference_length, insertions, code_args):
    logging.info("Getting node starts (%s starts, %s insertions", len(node_starts), len(insertions[0]))
    reference_node_lens = np.diff(np.r_[node_starts, reference_length])
    insertion_node_lens = insertions[1]
    node_lens = np.empty(reference_node_lens.size+insertion_node_lens.size)
    node_lens[code_args[-node_starts.size:]] = reference_node_lens
    node_lens[code_args[:-node_starts.size]] = insertion_node_lens
    return node_lens


def get_snps(snp_positions, reference_path, n_nodes, seqs=None, indices=None):
    logging.info("Creating SNPs %s", len(snp_positions))
    snp_positions = np.asanyarray(snp_positions, dtype="int")
    if not snp_positions.size:
        return SNPs()
    args = np.argsort(snp_positions)
    snp_positions = snp_positions[args]
    if seqs is not None:
        seqs = np.asanyarray(seqs)
        seqs = seqs[args]

    snp_positions.sort()
    snp_idxs = np.searchsorted(reference_path._distance_to_node, snp_positions, side="right")-1
    snp_offsets = snp_positions-reference_path._distance_to_node[snp_idxs]
    snp_node_ids = reference_path._node_ids[snp_idxs]
    node_index = np.searchsorted(snp_node_ids, np.arange(n_nodes))
    return SNPs(node_index, snp_offsets, seqs)


def build_seq_graph(insertion_seqs, insertion_node_ids, reference_path, fasta, n_nodes):
    logging.info("Building seq_graph")
    seqs = np.array([s.lower() for s in insertion_seqs])
    all_seqs = np.empty(n_nodes, dtype=object)
    all_seqs[insertion_node_ids] = seqs
    for node_id, start, stop in reference_path.get_node_intervals():
        all_seqs[node_id] = fasta[start:stop].seq.lower()
        assert len(all_seqs[node_id]) == stop-start, (start, stop, all_seqs[node_id])
    return all_seqs


def graph_from_snps_and_indels(deletions, insertions, snp_positions, reference_length, insertion_seqs=None, snp_seqs=None, fasta=None, vcf_map=None, insertion_var_ids=None):
    deletions = np.asanyarray(deletions, dtype="int")
    insertions = np.asanyarray(insertions, dtype="int")
    node_starts = get_node_starts(deletions, insertions)
    logging.info("Sorting")
    all_node_starts = np.concatenate((insertions[0], node_starts))
    tmp_code_args = np.argsort(all_node_starts, kind="mergesort")
    code_args = tmp_code_args.copy() # TODO: prob easyier way to to this
    code_args[tmp_code_args] = np.arange(tmp_code_args.size)
    reference_node_ids = code_args[-node_starts.size:]
    insertion_node_ids = code_args[:-node_starts.size]
    if vcf_map is not None:
        logging.info("Creating vcf map")
        vcf_map.fill_insertions(insertion_var_ids, insertion_node_ids)
    adj_list = create_adj_list(reference_node_ids, node_starts,
                               insertions, insertion_node_ids, deletions)
    node_lens = create_node_lens(node_starts, reference_length, insertions, code_args)
    reference_path = Path(reference_node_ids,
                          np.r_[0, np.cumsum(node_lens[reference_node_ids])])
    snps = get_snps(snp_positions, reference_path, all_node_starts.size, snp_seqs)
    if fasta is not None:
        seqs = build_seq_graph(insertion_seqs, insertion_node_ids, reference_path, fasta, all_node_starts.size)
    else:
        seqs = None
    return VCFGraph(node_lens, AdjList.from_dict(adj_list, all_node_starts.size), snps, seqs), reference_path


def construct_graphs(vcf_entries, reference_lengths, fasta):
    return ((chrom, construct_graph(entries, reference_lengths[int(chrom)],
                                    fasta[chrom]))
            for chrom, entries in vcf_entries if chrom.isdigit())


def construct_graph(vcf_entries, reference_length, fasta=None):
    if fasta is not None:
        assert len(fasta) == reference_length, (len(fasta), reference_length)
    insertion_positions = []
    insertion_lens = []
    deletion_starts = []
    deletion_lens = []
    snp_positions = []
    snp_seqs = []
    snp_var_ids = []
    ins_var_ids = []
    insertion_seqs = []
    i = 0
    counter = 0
    for entry in vcf_entries:
        if counter % 1000 == 0:
            logging.info("Entry %s" % counter)
        counter += 1
        if fasta:
            f = fasta[entry.pos:entry.pos+len(entry.ref)]
            assert str(f).lower() == entry.ref.lower(), (f, entry)
        if entry.alt.startswith("<"):
            continue
        var_type = classify_vcf_entry(entry)
        if var_type == SNP:
            new_entry = prune_SNP(entry)
            snp_positions.append(new_entry.pos)
            snp_seqs.append(new_entry.alt)
            snp_var_ids.append(i)
            i += 1
        elif var_type == INS:
            new_entry = prune_insertion(entry)
            insertion_positions.append(new_entry.pos)
            insertion_lens.append(len(new_entry.alt))
            insertion_seqs.append(new_entry.alt)
            ins_var_ids.append(i)
            i += 1
        elif var_type == DEL:
            new_entry = prune_deletion(entry)
            deletion_starts.append(new_entry.pos)
            deletion_lens.append(len(new_entry.ref))
    dels = np.array([deletion_starts, deletion_lens])
    insertions = np.array([insertion_positions, insertion_lens])
    logging.info("Finsihed parsing VCF")
    vcfmap = VCFMap.empty(i)
    vcfmap.fill_snps(snp_var_ids, snp_positions)
    graph, reference_path = graph_from_snps_and_indels(dels, insertions, snp_positions, reference_length, insertion_seqs, snp_seqs, fasta, vcf_map=vcfmap, insertion_var_ids=ins_var_ids)
    return graph, reference_path, vcfmap
