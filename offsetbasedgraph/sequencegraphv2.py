from .graph import BlockArray
import numpy as np
import json
import logging
import re


class SequenceGraphv2():
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

    def __init__(self, node_id_offset, sequence_array):
        self._node_id_offset = node_id_offset
        self._sequence_array = sequence_array
        self._sequence_cache = {}

    def to_file(self, file_name):
        file = open(file_name, "wb")
        np.savez_compressed(file,
                            node_id_offset=self._node_id_offset,
                            sequence_array=self._sequence_array)
        file.close()
    
    @classmethod
    def from_file(cls, file_name):
        file = open(file_name, "rb")
        data = np.load(file)
        seqgraph = cls(data["node_id_offset"],
                       data["sequence_array"])
        file.close()
        return seqgraph

    @classmethod
    def create_empty_from_ob_graph(cls, ob_graph):
        assert isinstance(ob_graph.blocks, BlockArray), "Block backend must be BlockArray"
        sequence_array = np.array(["X" * size for size in ob_graph.blocks._array], dtype=str)
        return cls(ob_graph.blocks.node_id_offset, sequence_array)

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

    def set_sequence(self, node_id, sequence):
        sequence = sequence.lower()
        new_sequence = re.sub(r'[^acgtn]', "n", sequence)
        if new_sequence != sequence:
            logging.warning("Trying to set invalid sequence to sequence graph: %s."
                            "Was replaced with n's, resulting in: %s" % (sequence, new_sequence))
            sequence = new_sequence

        index = node_id - self._node_id_offset
        try:
            pos = index
            self._sequence_array[pos] = sequence
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

    def get_node_sequence(self, node_id):
        index = node_id - self._node_id_offset
        pos = self._indices[index]
        node_size = self._node_sizes[index]
        sequence = self._sequence_array[pos:pos+node_size]
        return "".join(self._letters[sequence])

    def get_sequence(self, node_id, start=0, end=False):
        if start == 0 and not end:
            return self._sequence_array[node_id - self._node_id_offset]
        else:
            raise NotImplementedError("Not implemented, but should be trivial to do")


        if start == 0 and not end:
            cache_id = node_id
            if cache_id in self._sequence_cache:
                return self._sequence_cache[cache_id]

        index = node_id - self._node_id_offset
        sequence = self._sequence_array[index]

        assert len(sequence) == node_size

        if end is False:
            end = node_size
        sliced_sequence = sequence[start:end]
        if start == 0 and end == node_size:
            cache_id = node_id
            self._sequence_cache[cache_id] = sliced_sequence
        return sliced_sequence

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

