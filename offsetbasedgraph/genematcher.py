from collections import defaultdict
import pickle


class GeneMatcher(object):

    categories = ["FLANK", "VAR", "FLANK+VAR"]

    def __init__(self, alt_gene, main_genes):
        self.is_cut = False
        self.alt_gene = alt_gene
        self.category = self.classify_alt_gene(alt_gene)
        self.main_genes = [gene for gene in main_genes if
                           gene.strand == alt_gene.strand or
                           self.category == "VAR"]

        self.merged_rps = [rp for rp in
                           alt_gene.transcription_region.region_paths
                           if self.is_merged(rp)]
        self.find_match()

    def find_match(self):
        if not self.main_genes:
            self.scores = []
            self.score = -1
            return
        comp_func = self.compare_generic
        self.scores = [comp_func(main_gene) for main_gene in self.main_genes]
        max_score = max(self.scores)
        self.score = max_score
        self.match_gene = self.main_genes[
            self.scores.index(max_score)]

    @staticmethod
    def is_merged(name):
        return name.count("chr") > 1

    @staticmethod
    def was_alt(name):
        return "alt" in name

    @staticmethod
    def is_varying(name):
        return GeneMatcher.was_alt(name) and not GeneMatcher.is_merged(name)

    @staticmethod
    def classify_alt_gene(gene):
        rps = gene.transcription_region.region_paths
        if not any(GeneMatcher.is_varying(rp) for rp in rps):
            return "FLANK"
        if GeneMatcher.is_merged(rps[0]) and GeneMatcher.is_merged(rps[-1]):
            return "FLANK+VAR"
        elif GeneMatcher.is_merged(rps[0]):
            return "FLANK+VAR"
        elif GeneMatcher.is_merged(rps[-1]):
            return "FLANK+VAR"
        if len(rps) == 1:
            return "VAR"
        else:
            raise Exception("Couldnt classify %s" % rps)

    def compare_generic(self, main_gene):
        graph = self.alt_gene.graph
        alt_rps = self.alt_gene.transcription_region.region_paths
        main_rps = main_gene.transcription_region.region_paths
        paralell_rps = graph.find_parallell_blocks(alt_rps, graph.is_main_name)
        are_paralell = [rp in paralell_rps for rp in main_rps]
        d_length = abs(self.alt_gene.transcript_length -
                       main_gene.transcript_length)
        base_score = 0
        if all(are_paralell):
            if d_length < 5:
                return 5+base_score
            return 4+base_score
        if any(are_paralell) and main_gene.strand == self.alt_gene.strand:
            if d_length < 5:
                return 3+base_score
            return 2+base_score
        if d_length < 5:
            return 1

        return 0


class GeneMatchings(object):
    codes = {5: "Paralell and\t\t equal transcript length",
             4: "Paralell but\t\t unequal transcript length",
             3: "Partially paralell and\t equal transcript length",
             2: "Partially paralell but\t unequal transcript length",
             1: "Not paralell but\t equal transcript length",
             0: "Not paralell and\t unequal transcript length",
             -1: "No Matches on same strand\t\t \t"}

    def __init__(self, alt_genes, main_genes):
        self.alt_genes = alt_genes
        self.main_genes = main_genes
        self.matches = []
        self.find_matches()

    def find_matches(self):
        for name, alt_genes in self.alt_genes.lookup.items():
            main_genes = self.main_genes.lookup[name]
            for alt_gene in alt_genes:
                self.matches.append(GeneMatcher(alt_gene, main_genes))

    def calculate_summary_categories(self):
        cat_1 = len([m for m in self.matches
                     if m.category == "FLANK" and m.score == 5])
        cat_2 = len([m for m in self.matches
                     if m.category == "VAR" and m.score == 5])
        cat_3 = len([m for m in self.matches
                     if m.category == "FLANK+VAR" and m.score == 5])
        cat_4 = len([m for m in self.matches if m.score in (3, 2)])
        cat_5 = len(self.matches)-cat_1-cat_2-cat_3-cat_4
        return [cat_1, cat_2, cat_3, cat_4, cat_5]

    def find_unequal_length_uncut_isoforms(self):
        return [(m.alt_gene, m.match_gene) for m in
                self.matches if m.score == 2]

    def __str__(self):
        lines = []
        for category in GeneMatcher.categories:
            lines.append(category)
            category_mathes = [m for m in self.matches
                               if m.category == category]
            score_dict = defaultdict(int)
            for match in category_mathes:
                score_dict[match.score] += 1
            lines.extend("\t%s:\t %s" % (self.codes[k], v)
                         for k, v in score_dict.items())

        summary_categories = self.calculate_summary_categories()
        lines.extend(["", "Summary categories"])
        lines.extend(["Cat %s:\t %s" % (i+1, c)
                      for i, c in enumerate(summary_categories)])
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    def to_pickle(self, file_name):
        with open("%s" % file_name, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def from_pickle(cls, file_name):
        with open("%s" % file_name, "rb") as f:
            o = pickle.loads(f.read())
        assert isinstance(o, cls)
        return o
