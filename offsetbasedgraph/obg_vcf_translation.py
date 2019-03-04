from collections import namedtuple, defaultdict
import offsetbasedgraph as obg
from itertools import takewhile, chain
import logging
import pickle
import numpy as np

Pos = namedtuple("Position", ["node_id", "offset"])
SNP = namedtuple("SNP", ["snp_id"])
SNP_CODE = -1
Translation = namedtuple("Translation", ["node_id", "offset"])
VCFInterval = namedtuple("VCFInterval", ["start", "end", "node_ids", "snps"])


def var_type(variant):
    if len(variant.alt_seq) == 1 and len(variant.ref_seq) == 1:
        return "SNP"
    if len(variant.ref_seq) == 0:
        return "INS"
    if variant.alt_seq == variant.ref_seq:
        return "DUP"
    return None


class Translator:
    def __init__(self, translation, snp_index, extra_nodes, node_offset=0):
        self._translation = translation
        self._snp_index = snp_index
        self._extra_nodes = extra_nodes
        self._node_offset = node_offset

    def save(self, basename):
        np.save(basename+"_translation.npy", np.array(self._translation))
        np.save(basename+"_snp_index.npy", self._snp_index)
        np.save(basename+"_node_offset.npy", self._node_offset)
        pickle.dump(self._extra_nodes, open(basename+"_extra_nodes.pickle", "wb"))

    @classmethod
    def load(cls, basename):
        translation = Translation(*np.load(basename+"_translation.npy"))
        snp_index = np.load(basename+"_snp_index.npy")
        node_offset = np.load(basename+"_node_offset.npy")
        extra_nodes = pickle.load(open(basename+"_extra_nodes.pickle", "rb"))

        return cls(translation, snp_index, extra_nodes, node_offset)

    def translate_position(self, position, vcf_graph=None, is_start=True):
        node_id = self._translation.node_id[position.region_path_id]
        offset = self._translation.offset[position.region_path_id]+position.offset
        i = 0
        while (offset+is_start) > vcf_graph._node_lens[node_id]:
            offset -= vcf_graph._node_lens[node_id]
            node_id = self._extra_nodes[position.region_path_id][i]
            i += 1
        return Pos(node_id, offset)

    def get_snp_indices(self, node_ids):
        return [idx for idx in self._snp_index[node_ids] if idx >= 0 ]

    def offset_interval(self, interval):
        return obg.Interval(interval.start_position.offset,
                            interval.end_position.offset,
                            [rp-self._node_offset for rp in interval.region_paths],
                            graph=interval.graph)

    def _check_edges(self, node_ids, vcf_graph):
        for node_a, node_b in zip(node_ids[:-1], node_ids[1:]):
            if node_b not in vcf_graph._adj_list[node_a]:
                print("Node not in adj_list", node_a, node_b)

    def translate_interval(self, interval, vcf_graph):
        """TODO: Only handles forward intervals"""
        interval = self.offset_interval(interval)
        node_ids = self._translation.node_id[interval.region_paths]
        unique_node_ids = np.unique(node_ids)
        all_node_ids = []
        for i, node_id in enumerate(unique_node_ids):
            all_node_ids.append(node_id)
            if node_id in self._extra_nodes:
                all_node_ids.extend(self._extra_nodes[node_id])
        all_node_ids = np.unique(all_node_ids)
        self._check_edges(all_node_ids, vcf_graph)
        start = self.translate_position(
            interval.start_position, vcf_graph, True).offset
        end = self.translate_position(interval.end_position, vcf_graph, False)
        e_node = end.node_id
        end = end.offset
        while not e_node == all_node_ids[-1]:
            node_ids = all_node_ids[:-1]
        assert e_node == all_node_ids[-1], (e_node, all_node_ids)
        snps = self.get_snp_indices(interval.region_paths)
        new_interval = VCFInterval(start, end, all_node_ids, snps)
        return new_interval


class TranslationBuilder:
    def __init__(self, full_obg_graph, full_vcf_graph):
        self.full_obg_graph = full_obg_graph
        self.full_vcf_graph = full_vcf_graph
        size = full_obg_graph.graph._node_lens.size
        self.translation = Translation(
            -1*np.ones(size, dtype="int"),
            -1*np.ones(size, dtype="int"))
        self.snp_index = -1*np.ones(size, dtype="int")
        self.visited = np.zeros(size)
        self.extra_nodes = defaultdict(list)
        self._obg_sizes = self.full_obg_graph.graph._node_lens

    def build(self):
        logging.info("Building Ref tranlsation")
        self.build_ref_translation_numpy()
        logging.info("Building Alt translation")
        self.build_alt_translation()

        return Translator(self.translation, self.snp_index, self.extra_nodes, self.full_obg_graph.node_offset)

    def add_snp(self, variant):
        assert len(variant.ref_node_ids) == 1
        assert len(variant.alt_node_ids) == 1
        assert not any(self.full_obg_graph.path.is_in_path(node)
                       for node in variant.alt_node_ids)

        alt_node = variant.alt_node_ids[0]
        ref_node = variant.ref_node_ids[0]
        node_id = self.translation.node_id[ref_node]
        assert node_id != -1
        offset = self.translation.offset[ref_node]
        self.translation.node_id[alt_node] = node_id
        self.translation.offset[alt_node] = offset
        snp_id = self.full_vcf_graph.graph._snps.find_snp(node_id, offset, variant.alt_seq)
        if snp_id is not None:
            self.snp_index[variant.alt_node_ids[0]] = snp_id

    def add_ins(self, variant):
        assert not any(self.full_obg_graph.path.is_in_path(node)
                       for node in variant.alt_node_ids)

        node = variant.ref_node
        vcf_ref_node = self.translation.node_id[node]
        assert vcf_ref_node != -1
        offset = self.translation.offset[node]
        vcf_node_end = self.full_vcf_graph.path.distance_to_node_id(vcf_ref_node) + self.full_vcf_graph.graph._node_lens[vcf_ref_node]
        obg_node_end = self.full_obg_graph.path.distance_to_node_id(node) + self._obg_sizes[node]
        if not vcf_node_end == obg_node_end:
            print("Not registered INS at", vcf_node_end, obg_node_end)
            return False
        vcf_node = self.full_vcf_graph.find_insertion_from_node(
            vcf_ref_node, variant.alt_seq)
        assert vcf_node != -1
        self.add_translation(variant.alt_node_ids, vcf_node)
        return True

    def add_dup(self, variant):
        assert not any(self.full_obg_graph.path.is_in_path(node)
                       for node in variant.alt_node_ids)

        ref_node = variant.ref_node_ids[0]
        vcf_node = self.translation.node_id[ref_node]
        assert vcf_node != -1
        offset = self.translation.offset[ref_node]
        size = self.full_vcf_graph.graph._node_lens[vcf_node]
        for alt_node in variant.alt_node_ids:
            self.translation.node_id[alt_node] = vcf_node
            self.translation.offset[alt_node] = offset
            offset += self._obg_sizes[alt_node]
            while offset >= size:
                vcf_node = self.full_vcf_graph.next_node(vcf_node)
                if offset > size:
                    self.extra_nodes[alt_node].append(vcf_node)
                offset -= size
                size = self.full_vcf_graph.graph._node_lens[vcf_node]

    def build_alt_translation(self):
        obg_ref_nodes = self.full_obg_graph.traverse_ref_nodes()
        variants = chain.from_iterable(
            self.full_obg_graph.get_variants_from_node(int(node))
            for node in obg_ref_nodes)

        def handle_variant(variant):
            self.visited[variant.alt_node_ids] = 1
            t = var_type(variant)
            if t == "INS":
                self.add_ins(variant)
            elif t == "SNP":
                self.add_snp(variant)
            elif t == "DUP":
                self.add_dup(variant)
            else:
                assert False, (var_type(variant), variant)
        counter = 0
        for variant in variants:
            if counter % 1000 == 0:
                logging.info("Variant %s", counter)
            counter += 1
            handle_variant(variant)
        # assert (not res) or np.count_nonzero(self.translation.node_id[variant.alt_node_ids]==-1) == 0, (variant, t)

    def handle_multi_nodes(self, multi_nodes):
        print(multi_nodes)
        for m_node in multi_nodes:
            extra_nodes = []
            vcf_node_id = self.translation.node_id[m_node]
            offset = self.translation.offset[m_node]+self._obg_sizes[m_node]
            offset -= self.full_vcf_graph.graph._node_lens[vcf_node_id]
            while offset > 0:
                vcf_node_id = self.full_vcf_graph.next_node(vcf_node_id)
                extra_nodes.append(vcf_node_id)
                offset -= self.full_vcf_graph.graph._node_lens[vcf_node_id]
            self.extra_nodes[m_node] = extra_nodes

    def build_ref_translation_numpy(self):
        obg_distances = self.full_obg_graph.path._distance_to_node[:-1]
        obg_node_ids = self.full_obg_graph.path._node_ids
        vcf_distances = self.full_vcf_graph.path._distance_to_node
        obg_idxs = np.searchsorted(vcf_distances, obg_distances, side="right")-1
        vcf_node_idxs = self.full_vcf_graph.path._node_ids[obg_idxs]
        offsets = obg_distances-vcf_distances[obg_idxs]
        self.translation.node_id[obg_node_ids] = vcf_node_idxs
        self.translation.offset[obg_node_ids] = offsets

        multi_nodes = obg_distances + self._obg_sizes[obg_node_ids] > vcf_distances[obg_idxs+1]
        self.handle_multi_nodes(obg_node_ids[multi_nodes])

    def _build_ref_translation(self):
        obg_ref_nodes = self.full_obg_graph.traverse_ref_nodes()
        for vcf_node_id, start, end in self.full_vcf_graph.path.get_node_intervals():
            assert vcf_node_id != -1
            obg_nodes = []
            e = None
            for node in obg_ref_nodes:
                obg_nodes.append(node)
                s = self.full_obg_graph.path.distance_to_node_id(node)
                e = s+self._obg_sizes[node]
                if e == end:
                    break
                assert (e < end), (s, e, start, end, obg_nodes)
            else:
                print("FINISHED")
            assert e == end, (e, end)
            self.visited[obg_nodes] = 1
            obg_start = self.full_obg_graph.path.distance_to_node_id(obg_nodes[0])
            assert obg_start == start, (obg_start, start)
            obg_end = self.full_obg_graph.path.distance_to_node_id(
                obg_nodes[-1])+self._obg_sizes[obg_nodes[-1]]
            assert obg_end == e, (obg_end, end, obg_nodes)
            assert obg_end == end, (obg_end, end, obg_nodes)

            self.add_translation(obg_nodes, vcf_node_id)
            obg_seq = self.full_obg_graph.seq_graph.get_nodes_sequences(obg_nodes).lower()
            vcf_seq = self.full_vcf_graph.graph._seqs[vcf_node_id].lower()
            assert (obg_seq == vcf_seq), (obg_seq, vcf_seq)
            assert np.count_nonzero(self.translation.node_id[obg_nodes]==-1) == 0

    def build_ref_translation(self):
        obg_ref_nodes = self.full_obg_graph.traverse_ref_nodes()
        overlapping_node = None
        overlap_offset = 0
        for vcf_node_id, start, end in self.full_vcf_graph.path.get_node_intervals():
            assert vcf_node_id != -1
            if overlap_offset >= (end-start):
                self.extra_nodes[overlapping_node].append(vcf_node_id)
                overlap_offset = overlap_offset-(end-start)
                continue
            if overlapping_node is not None:
                self.extra_nodes[overlapping_node].append(vcf_node_id)
                # self.extra_nodes[overlapping_node].append(overlap_offset)
            obg_nodes = []
            e = None
            tmp_overlapping_node = None
            tmp_overlap_offset = 0
            for node in obg_ref_nodes:
                obg_nodes.append(node)
                s = self.full_obg_graph.path.distance_to_node_id(node)
                assert e is None or s == e, (s, e)
                e = s + self._obg_sizes[node]
                if e == end:
                    tmp_overlapping_node = None
                    tmp_overlap_offset = 0
                    break
                if e > end:
                    print("overflow", node, vcf_node_id, e, end)
                    tmp_overlapping_node = node
                    tmp_overlap_offset = int(e-end)
                    break
                assert (e < end), (s, e, start, end, obg_nodes)
            else:
                print("FINISHED")
            assert tmp_overlapping_node is not None or e == end, (e, end)
            self.visited[obg_nodes] = 1
            if overlapping_node is None:
                obg_start = self.full_obg_graph.path.distance_to_node_id(obg_nodes[0])
                assert obg_start == start, (obg_start, start)
                # obg_end = self.full_obg_graph.linear_path.get_offset_at_node(obg_nodes[-1])+self._obg_sizes[obg_nodes[-1]]
                # assert obg_end == e, (obg_end, end, obg_nodes)
                # assert obg_end == end, (obg_end, end, obg_nodes)

            self.add_translation(obg_nodes, vcf_node_id, overlap_offset)
            obg_seq = self.full_obg_graph.seq_graph.get_nodes_sequences(obg_nodes).lower()
            vcf_seq = self.full_vcf_graph.graph._seqs[vcf_node_id][overlap_offset:].lower()
            assert (obg_seq[:len(obg_seq)-tmp_overlap_offset] == vcf_seq), (obg_seq, vcf_seq, tmp_overlap_offset, self.full_vcf_graph.graph._node_lens[vcf_node_id])
            assert np.count_nonzero(self.translation.node_id[obg_nodes]==-1) == 0
            overlap_offset = tmp_overlap_offset
            overlapping_node = tmp_overlapping_node

    def add_translation(self, obg_nodes, vcf_node, start_size=0):
        sizes = [start_size] + [self._obg_sizes[node]
                                for node in obg_nodes[:-1]]
        self.translation.node_id[obg_nodes] = vcf_node
        offsets = np.cumsum(sizes)
        self.translation.offset[obg_nodes] = offsets
