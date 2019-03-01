from offsetbasedgraph.vcfmap import VCFEntry, INS, DEL, SNP, prune_insertion
import pytest


def test_start_prune():
    entry = VCFEntry(10, "cacac", "cacacag")
    pruned = prune_insertion(entry)
    assert pruned == VCFEntry(15, "", "ag")
