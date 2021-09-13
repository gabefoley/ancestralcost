import ancestral_cost.ancestral_cost as ac


# def test_get_positions_with_content_root(simple_tree_4, simple_aln_4):
#     positions = ac.get_positions_with_content(simple_tree_4, simple_aln_4, node="root")
#     assert positions == [0,2]

def test_get_positions_with_content_root(simple_tree_6, simple_aln_6):
    positions = ac.get_positions_with_content(simple_tree_6, simple_aln_6)
    assert positions == [0,2,3]

# def test_no_matching_label(simple_tree_6, simple_aln_6):
#     with pytest.raises(NameError):
#         positions = ac.get_positions_with_content(simple_tree_6, simple_aln_6, node_name="Won't find me")

def test_get_positions_with_content_labelled_position1(simple_tree_6, simple_aln_6):
    " Test re"
    positions = ac.get_positions_with_content(simple_tree_6, simple_aln_6, node_name="Labelled_node_1")
    assert positions == [0,1,2,3]


def test_get_positions_with_content_labelled_position2(simple_tree_6, simple_aln_6):
    " Test re"
    positions = ac.get_positions_with_content(simple_tree_6, simple_aln_6, node_name="Labelled_node_2")
    assert positions == [0,1,3]