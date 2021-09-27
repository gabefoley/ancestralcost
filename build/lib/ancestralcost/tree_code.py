from ete3 import PhyloTree
from Bio import AlignIO
import ancestralcost.indel_placement as indel_placement


def load_tree_with_alignment(tree, alignment):
    """Load a tree and associate it with a specific alignment"""
    tree = PhyloTree(tree, alignment=alignment, format=1, alg_format="fasta")
    return tree


def get_positions_with_content(tree, alignment, node_name="#N0"):
    """
    Return the alignment positions that must be present in the ancestor at the given node
    :param tree: The phylogenetic tree
    :param node: The node of interest (defaults to root node)
    :return:
    """

    selected_node = None

    # Get the node of interest
    if node_name in ("#N0", "N0"):
        selected_node = tree.get_tree_root()

    else:
        for node in tree.get_descendants():
            if node.name == node_name:
                selected_node = node
                break

        if selected_node is None:
            raise NameError("The node you selected is not in this tree")

    # Check each position in the alignment to see if any leaf in both children have content there
    # If a position appears in each child tree, no need to check it at higher nodes - it must appear at this node
    found_positions = []

    # If a position appears in neither child tree, no need to check it at higher nodes - it won't appear at this node
    skip_positions = []

    found_positions = sorted(
        get_found_positions(
            selected_node, len(alignment[0]), found_positions, skip_positions
        )
    )

    return found_positions


def get_found_positions(selected_node, aln_len, found_positions, skip_positions):
    """Just return the positions that are implied to be there"""

    # Get the leaf nodes under each child subtree
    child1_leaves = selected_node.children[0].get_leaves()
    child2_leaves = selected_node.children[1].get_leaves()

    for pos in range(aln_len):
        if pos not in found_positions and pos not in skip_positions:
            found1 = False
            found2 = False
            for leaf in child1_leaves:
                if leaf.sequence[pos] != "-":
                    found1 = True
                    break
            for leaf in child2_leaves:
                if leaf.sequence[pos] != "-":
                    found2 = True
                    break

            # Must appear (no need to search in higher nodes)
            if found1 and found2:
                found_positions.append(pos)

            # Doesn't appear (no need to search in higher nodes)
            elif not found1 and not found2:
                skip_positions.append(pos)

    if selected_node.is_root() or len(found_positions) + len(skip_positions) == aln_len:
        return found_positions

    else:
        return get_found_positions(
            selected_node.up, aln_len, found_positions, skip_positions
        )


def get_parsimony_table(tree):
    """Create the table that holds the parsimony scores"""
    print("Making total parsimony\n")
    parsimony = indel_placement.make_total_parsimony(tree)
    print("Sorting parsimony table\n")
    parsimony.sort_index(inplace=True)

    return parsimony


def get_parsimony_score(parsimony, node, pos, unbalanced=True):
    """Get a parsimony score for a specific position"""
    score = parsimony[pos][node]

    # Only return the score if cost for presence is higher than cost for absence
    if unbalanced:
        if score[0] > score[1]:
            return score

    else:
        return score


def get_parsimony_scores(parsimony, positions, node="#N0"):
    """Get parsimony scores for all positions in a given node"""
    parsimony_dict = {}
    for pos in positions:
        score = get_parsimony_score(parsimony, node, pos)
        if score:
            parsimony_dict[pos] = score
    return parsimony_dict
