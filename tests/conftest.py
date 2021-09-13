import pytest
from Bio import AlignIO


import ancestral_cost.ancestral_cost as ac
file_folder = "./files/"


@pytest.fixture()
def simple_tree_4():

    simple_tree_4 = ac.load_tree_with_alignment(file_folder + "simple_4.nwk", file_folder + "simple_4.aln")
    return simple_tree_4

@pytest.fixture()
def simple_tree_6():
    simple_tree_6 = ac.load_tree_with_alignment(file_folder + "simple_6.nwk", file_folder + "simple_6.aln")
    return simple_tree_6

@pytest.fixture()
def simple_aln_4():
    simple_aln_4 = AlignIO.read(file_folder + "simple_4.aln", "fasta")
    return simple_aln_4

@pytest.fixture()
def simple_aln_6():
    simple_aln_6 = AlignIO.read(file_folder + "simple_6.aln", "fasta")
    return simple_aln_6
