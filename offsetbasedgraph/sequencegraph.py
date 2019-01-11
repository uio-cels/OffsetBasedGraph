from .graph import BlockArray
import numpy as np
import json
import logging
import re


class SequenceGraph():
    _compliments = {"A": "T",
                    "a": "t",
                    "C": "G",
                    "c": "g",
                    "T": "A",
                    "t": "a",
                    "G": "C",
                    "g": "c",
                    "n": "n",
                    "N": "N"}

    _letters = np.array(["n", "a", "c", "t", "g", "m"])

    def __init__(self, node_id_offset, indices, node_sizes, sequence_array):
        self._indices = indices
        self._node_id_offset = node_id_offset
        self._node_sizes = node_sizes
        self._sequence_array = sequence_array

    def to_file(self, file_name):
        file = open(file_name, "wb")
        np.savez_compressed(file,
                            node_id_offset=self._node_id_offset,
                            indices=self._indices,
                            node_sizes=self._node_sizes,
                            sequence_array=self._sequence_array)
        file.close()
    
    @classmethod
    def from_file(cls, file_name):
        file = open(file_name, "rb")
        data = np.load(file)
        seqgraph = cls(data["node_id_offset"],
                       data["indices"],
                       data["node_sizes"],
                       data["sequence_array"])
        file.close()
        return seqgraph

    @classmethod
    def create_empty_from_ob_graph(cls, ob_graph):
        assert isinstance(ob_graph.blocks, BlockArray), "Block backend must be BlockArray"
        blocks = ob_graph.blocks
        node_sizes = ob_graph.blocks._array
        indices = np.cumsum(node_sizes)
        indices = np.insert(indices, 0, 0)
        sequence_array = np.zeros(np.sum(node_sizes), dtype=np.uint8)
        return cls(ob_graph.blocks.node_id_offset, indices, node_sizes, sequence_array)

    def set_sequences_using_vg_json_graph(self, file_name):
        logging.info("Setting sequences using vg json graph %s" % file_name)
        f = open(file_name)
        i = 1
        for line in f.readlines():
            line = json.loads(line)
            if "node" in line:
                node_objects = line["node"]
                for node_object in node_objects:
                    i += 1
                    if i % 1000000 == 0:
                        logging.info("%d nodes processed" % i)
                    self.set_sequence(int(node_object["id"]), node_object["sequence"])

    def _letter_sequence_to_numeric(self, sequence):
        assert isinstance(sequence, np.ndarray), "Sequence must be numpy array"
        numeric = np.zeros_like(sequence, dtype=np.uint8)
        numeric[np.where(sequence == "n")[0]] = 0
        numeric[np.where(sequence == "a")[0]] = 1
        numeric[np.where(sequence == "c")[0]] = 2
        numeric[np.where(sequence == "t")[0]] = 3
        numeric[np.where(sequence == "g")[0]] = 4
        numeric[np.where(sequence == "m")[0]] = 4
        return numeric

    def set_sequence(self, node_id, sequence):
        sequence = sequence.lower()
        new_sequence = re.sub(r'[^acgtn]', "n", sequence)
        if new_sequence != sequence:
            logging.warning("Trying to set invalid sequence to sequence graph: %s."
                            "Was replaced with n's, resulting in: %s" % (sequence, new_sequence))
            sequence = new_sequence

        index = node_id - self._node_id_offset
        node_size = self._node_sizes[index]
        assert node_size == len(sequence), "Invalid sequence. Does not cover whole node"
        sequence = np.array(list(sequence))
        sequence = self._letter_sequence_to_numeric(sequence)
        try:
            pos = self._indices[index]
            self._sequence_array[pos:pos+node_size] = sequence
        except ValueError:
            raise ValueError("Trying to create sequence graph with "
                             "invalid sequence (not containing A,C,G,T,N): %s" % sequence)

    def get_sequence_on_directed_node(self, node_id, start=0, end=False):
        """Handles directed nodes"""
        if node_id > 0:
            return self.get_sequence(node_id, start, end)
        else:
            node_size = len(self.get_sequence(-node_id))

            if not end:
                end = node_size

            new_start = node_size - end
            new_end = node_size - start
            reverse_sequence = self.get_sequence(-node_id, new_start, new_end)
            return self._reverse_compliment(reverse_sequence)

    def get_numeric_node_sequence(self, node_id):
        index = node_id - self._node_id_offset
        pos = self._indices[index]
        node_size = self._node_sizes[index]
        sequence = self._sequence_array[pos:pos+node_size]
        return sequence

    def get_node_sequence(self, node_id):
        index = node_id - self._node_id_offset
        pos = self._indices[index]
        node_size = self._node_sizes[index]
        sequence = self._sequence_array[pos:pos+node_size]
        return "".join(self._letters[sequence])

    def get_sequence(self, node_id, start=0, end=False):
        index = node_id - self._node_id_offset
        pos = self._indices[index]
        node_size = self._node_sizes[index]
        sequence = self._sequence_array[pos:pos+node_size]

        assert len(sequence) == node_size

        if end is False:
            end = node_size
        return "".join(self._letters[sequence[start:end]])
        assert end <= node_size, "End is > node size (%d, %d)" % (end, node_size)

        assert end >= start, "End in get_sequence is smaller than start (start: %d, end: %d. Node size: %d). Node: %d. Node id offset: %d" % (start, end, node_size, node_id, self._node_id_offset)

        subsequence = ''.join(self._letters[sequence[start:end]])
        assert len(subsequence) == end - start, "len subsequence: %d, end:%d, start%d, node: %d. Node size: %d" % (len(subsequence), end, start, node_id, node_size)
        return subsequence

    def _reverse_compliment(self, sequenence):
        return "".join(self._compliments[c] for c in sequenence[::-1])

    def node_size(self, node_id):
        return self._node_sizes[abs(node_id)-self._node_id_offset]

    def get_interval_sequence(self, interval):

        rps = interval.region_paths
        start_node = rps[0]
        end_node = rps[-1]

        if start_node == end_node and len(rps) == 1:
            ret_str =  self.get_sequence_on_directed_node(
                start_node,
                interval.start_position.offset,
                interval.end_position.offset)
        else:

            start_sequence = self.get_sequence_on_directed_node(
                start_node,
                interval.start_position.offset)
            end_sequence = self.get_sequence_on_directed_node(
                end_node,
                0,
                interval.end_position.offset)
            middle_sequence = ""
            for rp in rps[1:-1]:
                middle_sequence += self.get_sequence_on_directed_node(rp)

            ret_str =  "%s%s%s" % (start_sequence, middle_sequence, end_sequence)
        if interval.graph is None:
            interval.graph = self
        assert len(ret_str) == interval.length(), (ret_str, interval)
        return ret_str

    def count_gc_content(self):
        alphabet = np.array([[1], [2], [3], [4]])
        counts = np.count_nonzero(self._sequence_array == alphabet, axis=1)
        g_c = counts[1]+counts[3]
        a_t = counts[0]+counts[2]
        return self._sequence_array.size, g_c, a_t
