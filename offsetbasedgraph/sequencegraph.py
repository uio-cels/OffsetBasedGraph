from .graph import BlockArray
import numpy as np
import json
import logging

class SequenceGraph():
    _compliments = {"A": "T",
                    "a": "t",
                    "C": "G",
                    "c": "g",
                    "T": "A",
                    "t": "a",
                    "G": "C",
                    "g": "c"}

    _letters = np.array(["n", "a", "c", "t", "g"])

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
        sequence[np.where(sequence == "n")[0]] = 0
        sequence[np.where(sequence == "a")[0]] = 1
        sequence[np.where(sequence == "c")[0]] = 2
        sequence[np.where(sequence == "t")[0]] = 3
        sequence[np.where(sequence == "g")[0]] = 4
        return sequence

    def set_sequence(self, node_id, sequence):
        sequence = sequence.lower()
        index = node_id - self._node_id_offset
        node_size = self._node_sizes[index]
        assert node_size == len(sequence), "Invalid sequence. Does not cover whole node"
        sequence = np.array(list(sequence))
        sequence = self._letter_sequence_to_numeric(sequence)
        pos = self._indices[index]
        self._sequence_array[pos:pos+node_size] = sequence

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

    def get_sequence(self, node_id, start=0, end=False):
        index = node_id - self._node_id_offset
        pos = self._indices[index]
        node_size = self._node_sizes[index]
        sequence = self._sequence_array[pos:pos+node_size]

        if not end:
            end = node_size

        return ''.join(self._letters[sequence[start:end]])

    def _reverse_compliment(self, sequenence):
        return "".join(self._compliments[c] for c in sequenence[::-1])

    def get_interval_sequence(self, interval):

        rps = interval.region_paths
        start_node = rps[0]
        end_node = rps[-1]

        if start_node == end_node and len(rps) == 1:
            return self.get_sequence_on_directed_node(
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

            return "%s%s%s" % (start_sequence, middle_sequence, end_sequence)