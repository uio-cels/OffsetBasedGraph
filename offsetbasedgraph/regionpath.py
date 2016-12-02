class RegionPath(object):
    """Region Path as described in the article.
    Contains reference to the areas of the linear
    reference genome that it covers
    """

    def __init__(self, _id, linear_reference_dict):
        self.id = _id
        self.linear_references = linear_reference_dict
        self.main_path_linear_reference = None

    def get_length(self):
        return max([lr.end-lr.start for lr in self.linear_references.values()])

    def contains(self, lin_ref):
        return any([my_ref.contains(lin_ref) for my_ref
                    in self.linear_references.values()])

    def distance_to_position(self, chrom, offset):
        assert chrom in self.linear_references
        lin_ref = self.linear_references[chrom]
        return lin_ref.distance_to_position(chrom, offset)

    def contains_position(self, chrom, offset):
        if chrom not in self.linear_references:
            return False
        return self.linear_references[chrom].contains_position(chrom, offset)

    def intersects(self, lin_ref):
        return any([my_ref.intersects(lin_ref) for my_ref
                    in self.linear_references.values()])

    def length(self):
        return max([lr.length() for lr in self.linear_references.values()])

    def is_empty(self):
        if not self.linear_references:
            return True
        if (all([lr.start == lr.end for lr
                 in self.linear_references.values()])):
            return True
        return False

    def split(self, offset):
        """Split the region path at the given offsets"""
        assert offset < self.length()
        lin_ref_dicts = {key: val.split() for
                         key, val in self.linear_references.items()}

        lin_ref_dicts = [{key: val[i] for key, val in lin_ref_dicts.items()}
                         for i in range(len(offset))]

        ids = [self.id] + [self.id+str(i) for i in range(1, len(offset))]

        return [RegionPath(ids[i], lin_ref_dict)
                for i, lin_ref_dict in enumeratei(lin_ref_dicts)]

    def __eq__(self, other):
        self_refs = self.linear_references
        other_refs = other.linear_references

        if not len(self_refs) == len(other_refs):
            return False

        for lr in self_refs:
            if lr not in other_refs:
                return False

        return True

    def __str__(self):
        return " RegionPath: %s"%self.id + ";".join([str(lin_ref) for lin_ref
                                                     in self.linear_references.items()])

    def __repr__(self):
        return self.__str__()
