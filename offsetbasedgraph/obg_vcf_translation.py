from collections import namedtuple, defaultdict
from itertools import takewhile, chain
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
    def __init__(self, translation, snp_index, extra_nodes):
        self._translation = translation
        self._snp_index = snp_index
        self._extra_nodes = extra_nodes

    def save(self, basename):
        np.save(basename+"_translation.npy", np.array(self._translation))
        np.save(basename+"_snp_index.npy", self._snp_index)
        pickle.dump(self._extra_nodes, open(basename+"_extra_nodes.pickle", "wb"))

    @classmethod
    def load(cls, basename):
        translation = Translation(*np.load(basename+"_translation.npy"))
        snp_index = np.load(basename+"_snp_index.npy")
        extra_nodes = pickle.load(open(basename+"_extra_nodes.pickle", "rb"))
        return cls(translation, snp_index, extra_nodes)

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

    def translate_interval(self, interval, vcf_graph):
        """TODO: Only handles forward intervals"""

        node_ids = self._translation.node_id[interval.region_paths]
        offsets = self._translation.offset[interval.region_paths]
        changes = np.flatnonzero(np.diff(node_ids))
        node_sizes = np.array([interval.graph.node_size(n) for n in interval.region_paths])

        unique_node_ids = np.unique(node_ids)
        ends = offsets[changes]+node_sizes[changes]
        e_node_ids = []
        for i, node_id in enumerate(unique_node_ids):
            e_node_ids.append(node_id)
            if node_id in self._extra_nodes:
                e_node_ids.extend(self._extra_nodes[node_id])
        unique_node_ids = np.array(e_node_ids, dtype="int")
        for node_a, node_b in zip(unique_node_ids[:-1], unique_node_ids[1:]):
            if node_b not in vcf_graph._adj_list[node_a]:
                print("Node not in adj_list")
                print(node_a, node_b, vcf_graph._adj_list[node_a], interval)
                print(node_ids, unique_node_ids)
        # true_ends = vcf_graph._node_lens[unique_node_ids[:-1]]
        # assert np.all(ends == true_ends), (ends, true_ends, node_ids, offsets, interval)
        # assert np.all(offsets[changes+1] == 0), (node_ids, offsets, interval)
        node_ids = np.unique(node_ids)

        # print(node_ids)
        start = self.translate_position(interval.start_position, vcf_graph, True).offset
        end = self.translate_position(interval.end_position, vcf_graph, False)
        e_node = end.node_id
        end = end.offset 
        if not e_node == node_ids[-1]:
            node_ids = node_ids[:-1]
        assert e_node == node_ids[-1], (e_node, node_ids)

        snps = self.get_snp_indices(interval.region_paths)
        new_interval = VCFInterval(start, end, node_ids, snps)
        # assert vcf_graph.interval_length(new_interval) == interval.length(), (interval, new_interval)
        return new_interval


class TranslationBuilder:
    def __init__(self, full_obg_graph, full_vcf_graph):
        self.full_obg_graph = full_obg_graph
        self.full_vcf_graph = full_vcf_graph
        size = full_obg_graph.graph.blocks.max_node_id()+1
        self.translation = Translation(
            -1*np.ones(size, dtype="int"),
            -1*np.ones(size, dtype="int"))
        self.snp_index = -1*np.ones(size, dtype="int")
        self.visited = np.zeros(size)
        self.extra_nodes = defaultdict(list)

    def build(self):
        self.build_ref_translation()
        for node_id in self.full_obg_graph.traverse_ref_nodes():
            obg_pos = self.full_obg_graph.linear_path.get_offset_at_node(node_id)
            vcf_node = self.translation.node_id[node_id]
            vcf_pos = self.full_vcf_graph.path.distance_to_node_id(vcf_node)
            assert vcf_pos<=obg_pos<=vcf_pos+self.full_vcf_graph.graph._node_lens[vcf_node], (node_id, obg_pos, vcf_node, vcf_pos)
        print("EXTRA:")
        print(self.extra_nodes)
        print("------------")
        assert not any(self.translation.node_id[node_id]== -1 for node_id in self.full_obg_graph.traverse_ref_nodes())
        self.build_alt_translation()
        for node_id in self.full_obg_graph.traverse_ref_nodes():
            obg_pos = self.full_obg_graph.linear_path.get_offset_at_node(node_id)
            vcf_node = self.translation.node_id[node_id]
            vcf_pos = self.full_vcf_graph.path.distance_to_node_id(vcf_node)
            assert vcf_pos<=obg_pos<=vcf_pos+self.full_vcf_graph.graph._node_lens[vcf_node], (node_id, obg_pos, vcf_node, vcf_pos)

        return Translator(self.translation, self.snp_index, self.extra_nodes)

    def add_snp(self, variant):
        assert len(variant.ref_node_ids) == 1
        assert len(variant.alt_node_ids) == 1
        assert not any(node in self.full_obg_graph.linear_path.nodes_in_interval()
                       for node in variant.alt_node_ids)
        assert not all(node in self.full_obg_graph.linear_path.nodes_in_interval()
                       for node in variant.alt_node_ids)

        alt_node = variant.alt_node_ids[0]
        ref_node = variant.ref_node_ids[0]
        node_id = self.translation.node_id[ref_node]
        assert node_id != -1
        offset = self.translation.offset[ref_node]
        self.translation.node_id[alt_node] = node_id
        if node_id == 14081:
            print("NODE", 14081)
            print(alt_node, node_id, ref_node)
        self.translation.offset[alt_node] = offset
        snp_id = self.full_vcf_graph.graph._snps.find_snp(node_id, offset, variant.alt_seq)
        self.snp_index[variant.alt_node_ids[0]] = snp_id

    def add_ins(self, variant):
        assert not any(node in self.full_obg_graph.linear_path.nodes_in_interval()
                       for node in variant.alt_node_ids)
        assert not all(node in self.full_obg_graph.linear_path.nodes_in_interval()
                       for node in variant.alt_node_ids)

        node = variant.ref_node
        vcf_ref_node = self.translation.node_id[node]
        offset = self.translation.offset[node]
        vcf_node_end = self.full_vcf_graph.path.distance_to_node_id(vcf_ref_node) + self.full_vcf_graph.graph._node_lens[vcf_ref_node]
        obg_node_end = self.full_obg_graph.linear_path.get_offset_at_node(node) + self.full_obg_graph.graph.node_size(node)
        if not vcf_node_end == obg_node_end:
            print("Not registered INS at", vcf_node_end, obg_node_end)
            return False
        vcf_node = self.full_vcf_graph.find_insertion_from_node(
            vcf_ref_node, variant.alt_seq)
        assert vcf_node != -1
        self.add_translation(variant.alt_node_ids, vcf_node)
        return True

    def add_dup(self, variant):
        assert not any(node in self.full_obg_graph.linear_path.nodes_in_interval()
                       for node in variant.alt_node_ids)
        assert not all(node in self.full_obg_graph.linear_path.nodes_in_interval()
                       for node in variant.alt_node_ids)

        ref_node = variant.ref_node_ids[0]
        vcf_node = self.translation.node_id[ref_node]
        assert vcf_node != -1
        offset = self.translation.offset[ref_node]
        size = self.full_vcf_graph.graph._node_lens[vcf_node]
        for alt_node in variant.alt_node_ids:
            self.translation.node_id[alt_node] = vcf_node
            self.translation.offset[alt_node] = offset
            offset += self.full_obg_graph.graph.node_size(alt_node)
            while offset >= size:
                vcf_node = self.full_vcf_graph.next_node(vcf_node)
                if offset > size:
                    self.extra_nodes[alt_node].append(vcf_node)
                offset -= size
                size = self.full_vcf_graph.graph._node_lens[vcf_node]

    def build_alt_translation(self):
        obg_ref_nodes = self.full_obg_graph.traverse_ref_nodes()
        for node in obg_ref_nodes:
            node = int(node)
            for variant in self.full_obg_graph.get_variants_from_node(node):
                res = True
                self.visited[variant.alt_node_ids] = 1
                t = var_type(variant)
                if t == "INS":
                    res = self.add_ins(variant)
                elif t == "SNP":
                    self.add_snp(variant)
                elif t == "DUP":
                    self.add_dup(variant)
                else:
                    assert False, (var_type, variant)
                assert (not res) or np.count_nonzero(self.translation.node_id[variant.alt_node_ids]==-1) == 0, (variant, t)

    def build_ref_translation_numpy(self):
        obg_distances = self.full_obg_graph.linear_path._node_to_distance
        idxs = np.flatnonzero(obg_distances[:-1] != obg_distances[1:])+1
        obg_distances = obg_distances[idxs]
        vcf_distances = self.full_vcf_graph.path._distance_to_node
        obg_idxs = np.searchsorted(vcf_distances, obg_distances, side="right")-1
        vcf_node_idxs = self.full_vcf_graph.path._node_ids[obg_idxs]
        offsets = obg_distances-vcf_distances[obg_idxs]
        self.translation.node_id[idxs+1] = vcf_node_idxs
        self.translation.node_id[idxs+1] = offsets

        multi_nodes = obg_distances > vcf_distances[obg_idxs+1]
        self.handle_multi_nodes(idxs[multi_nodes]+1)

    def _build_ref_translation(self):
        obg_ref_nodes = self.full_obg_graph.traverse_ref_nodes()
        for vcf_node_id, start, end in self.full_vcf_graph.path.get_node_intervals():
            assert vcf_node_id != -1
            # print(start, end)
            obg_nodes = []
            e = None
            for node in obg_ref_nodes:
                obg_nodes.append(node)
                s = self.full_obg_graph.linear_path.get_offset_at_node(node)
                e = s+self.full_obg_graph.graph.node_size(node)
                if e == end:
                    break
                assert (e < end), (s, e, start, end, obg_nodes)
            else:
                print("FINISHED")
            assert e == end, (e, end)
            self.visited[obg_nodes] = 1
            # obg_nodes = list(takewhile(lambda node: self.full_obg_graph.linear_path.get_offset_at_node(node) < end, obg_ref_nodes))
            obg_start = self.full_obg_graph.linear_path.get_offset_at_node(obg_nodes[0])
            assert obg_start == start, (obg_start, start)
            obg_end = self.full_obg_graph.linear_path.get_offset_at_node(obg_nodes[-1])+self.full_obg_graph.graph.node_size(obg_nodes[-1])
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
                s = self.full_obg_graph.linear_path.get_offset_at_node(node)
                assert e is None or s == e, (s, e)
                e = s + self.full_obg_graph.graph.node_size(node)
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
                obg_start = self.full_obg_graph.linear_path.get_offset_at_node(obg_nodes[0])
                assert obg_start == start, (obg_start, start)
                # obg_end = self.full_obg_graph.linear_path.get_offset_at_node(obg_nodes[-1])+self.full_obg_graph.graph.node_size(obg_nodes[-1])
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
        sizes = [start_size] + [self.full_obg_graph.graph.node_size(node)
                                for node in obg_nodes[:-1]]
        self.translation.node_id[obg_nodes] = vcf_node
        offsets = np.cumsum(sizes)
        self.translation.offset[obg_nodes] = offsets


# class Translation:
#     def __init__(self, node_ids, offsets):
#         self._node_ids = node_ids
#         self._offsets = offsets
# 
#     def translate_node_id(self, node_id):
#         vcf_node_id = self._node_ids[node_id]
#         if vcf_node_id == SNP_CODE:
#             return SNP(self._offsets[node_id])
#         return Pos(vcf_node_id, self._offsets[node_id])
# 
#     def translate_position(self, position):
#         pos = self.translate_node_id(position.node_id)
#         if isinstance(pos, SNP):
#             return pos
#         return Pos(pos.node_id, pos.offset+position.offset)
# 
#     @classmethod
#     def add_link(self, obg_nodes, vcf_node, translation):
#         translation.node_id[obg_nodes] = vcf_node
#         offsets = list(accumulate(chain([0], (
#         translation.offset
#         
#         
# 
#     @classmethod
#     def create_from_graphs(cls, full_obg_graph, full_vcf_graph):
#         translation = np.zeros((full_obg_graph.blocks.n_nodes, 2))
#         fill_ref_translation(full_obg_graph, full_vcf_graph, translation)
#         for ref_node in full_obg_graph.traverse_ref_nodes():
#             variants = full_obg_graph.get_variants_from_node(ref_node)
#             vcf_ref_node, vcf_offset = translation[ref_node]
#             for variant in variants:
#                 t = var_type(variant)
#                 if t == "INS":
#                     vcf_node = full_vcf_graph.find_insertion_from_node(vcf_ref_node, variant.alt_seq)
#                     translation[variant.alt_nodes, 0] = vcf_node
#                     translation[variant.alt_nodes] = np.cumsum(
# 
#     @classmethod
#     def create_from_maps(obg_to_vcf_map, vcf_to_graph_map, translation):
#         insertions = obg_to_vcf_map.ins_map != 0
#         translation[insertions, 0] = vcf_to_graph_map._ids[obg_to_vcf_map.ins_map[insertions]]
#         translation[insertions, 1] = vcf_to_graph_map._offsets[obg_to_vcf_map.ins_map[insertions]]
#         snps = obg_to_vcf_map.snp_map != 0
#         translation[snps, 0] = SNP_CODE
#         translation[snps, 1] = vcf_to_graph_map._ids[obg_to_vcf_map[snps]]
# 
# def fill_ref_translation(full_graph, full_vcf_graph, translation):
#     lookup = {}
#     cur_vcf_node_idx = 0
#     cur_vcf_offset = 0
#     path = full_vcf_graph.path
#     vcf_graph = full_vcf_graph.graph
#     for obg_node_id in full_graph.linear_path.get_sorted_nodes_in_interval():
#         offset = int(full_graph.linear_path.get_offset_at_node(obg_node_id))
#         if offset >= path._distance_to_node[cur_vcf_node_idx+1]:
#             cur_vcf_node_idx += 1
#             cur_vcf_offset = path._distance_to_node[cur_vcf_node_idx]
#             assert cur_vcf_offset == offset, (cur_vcf_offset, offset)
#         vcf_node_id = path._node_ids[cur_vcf_node_idx]
#         translation[obg_node_id] = (vcf_node_id, offset-cur_vcf_offset)
#         obg_seq = full_graph.seq_graph.get_sequence_on_directed_node(obg_node_id)
#         vcf_node_offset = offset-cur_vcf_offset
#         vcf_seq = vcf_graph._seqs[vcf_node_id][vcf_node_offset:vcf_node_offset+len(obg_seq)]
#         assert obg_seq == vcf_seq, (obg_seq, vcf_seq, obg_node_id, vcf_node_id)
# 
# def get_translation(ob_graph, sequence_graph, linear_interval,  vcf_graph, path):
#     lookup = {}
#     cur_vcf_node_idx = 0
#     cur_vcf_offset = 0
#     for obg_node_id in linear_interval.get_sorted_nodes_in_interval():
#         offset = int(linear_interval.get_offset_at_node(obg_node_id))
#         if offset >= path._distance_to_node[cur_vcf_node_idx+1]:
#             cur_vcf_node_idx += 1
#             cur_vcf_offset = path._distance_to_node[cur_vcf_node_idx]
#             assert cur_vcf_offset == offset, (cur_vcf_offset, offset)
#         vcf_node_id = path._node_ids[cur_vcf_node_idx]
#         lookup[obg_node_id] = (vcf_node_id, offset-cur_vcf_offset)
#         obg_seq = sequence_graph.get_sequence_on_directed_node(obg_node_id)
#         vcf_node_offset = offset-cur_vcf_offset
#         print("#", vcf_node_offset)
#         vcf_seq = vcf_graph._seqs[vcf_node_id][vcf_node_offset:vcf_node_offset+len(obg_seq)]
#         assert obg_seq == vcf_seq, (obg_seq, vcf_seq, obg_node_id, vcf_node_id)
#     return lookup
# 
# 
