from .vcfmap import *
from collections import defaultdict
import numpy as np
import logging

char_codes = {"a": 0, "c": 1, "g": 2, "t": 3}
char_rev = np.array(["a", "c", "g", "t"])


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
    def __init__(self, node_index=[0], snps=[], seqs=None, char_seqs=True):
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


class AdjList:
    def __init__(self, node_index=[0], to_nodes=[]):
        self._node_index = np.asanyarray(node_index, dtype="int")
        self._to_nodes = np.asanyarray(to_nodes, dtype="int")

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

    def __repr__(self):
        return "AdjList(%s, %s)" % (self._node_index, self._to_nodes)


class VCFGraph:
    def __init__(self, node_lens, adj_list, snps, seqs=None):
        self._node_lens = np.asanyarray(node_lens, dtype="int")
        self._adj_list = adj_list
        self._snps = snps
        self._seqs = np.asanyarray(seqs)

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
    def from_vcf(filename):
        with open(filename) as f:
            pass


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
    seqs = np.array(insertion_seqs)
    all_seqs = np.empty(n_nodes, dtype=object)
    all_seqs[insertion_node_ids] = seqs
    for node_id, start, stop in reference_path.get_node_intervals():
        all_seqs[node_id] = fasta[start:stop].seq
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
