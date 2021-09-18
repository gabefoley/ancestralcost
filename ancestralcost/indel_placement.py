from collections import defaultdict
import pandas as pd
from copy import deepcopy
import ancestralcost.alignment_permutations as alignment_permutations
from ete3 import PhyloTree, TreeStyle, TextFace, add_face_to_node, SeqMotifFace
from Bio import AlignIO


def fill_in_phylo_contigs(phylo_contigs, length, skip=[]):
    complete_phylo_contigs = defaultdict(list)
    # If there are no restrictions on which positions we have to add into (i.e. the sequence is blank at this point,
    # we can add in any combination

    if phylo_contigs == []:
        complete_phylo_contigs["rmv"].append([x for x in range(0, length)])
        complete_phylo_contigs["add"].append([x for x in range(0, length)])

    for idx, pos in enumerate(phylo_contigs):
        present = []
        if idx + 1 <= len(phylo_contigs):
            if pos != 0:
                complete_phylo_contigs["rmv"].append(
                    [x for x in range(0, pos) if x not in skip]
                )

            present.append(pos)
            while idx + 1 < len(phylo_contigs):
                next = phylo_contigs[idx + 1]

                if next == pos + 1:
                    present.append(next)
                else:
                    complete_phylo_contigs["add"].append(
                        [x for x in present if x not in skip]
                    )
                    complete_phylo_contigs["rmv"].append(
                        [x for x in range(present[-1] + 1, next) if x not in skip]
                    )
                    present = []

                pos = next
                if next not in present:
                    present.append(next)
                idx += 1
            if len(present) > 0:
                complete_phylo_contigs["add"].append(
                    [x for x in present if x not in skip]
                )
            if pos != length and pos != length - 1:
                complete_phylo_contigs["rmv"].append(
                    [x for x in range(pos + 1, length) if x not in skip]
                )
            break

    if len(complete_phylo_contigs) == 0:
        complete_phylo_contigs["rmv"].append([x for x in range(0, length)])
        complete_phylo_contigs["add"].append([x for x in range(0, length)])

    return complete_phylo_contigs


# fill_in_phylo_contigs([0,1,3], 5)

# fill_in_phylo_contigs([0,1,2,3,7,8,9], 10)
#
fill_in_phylo_contigs([2, 3, 6, 7, 12, 14, 19], 22)


def label_internal_nodes(tree):
    """
    Add internal labels to a tree. Returns a dictionary of previous labels mapping from previous label to new label
    """
    previous_labels = {}
    idx = 0
    for node in tree.traverse("preorder"):
        if not node.is_leaf():
            if node.name:
                previous_labels[node.name] = "#N" + str(idx)
            node.name = "#N" + str(idx)
            idx += 1

    return previous_labels


def make_extant_parsimony(tree):
    """
    Make a parsimony table for the extant sequences
    :param tree: Tree with annotated sequences
    :return: A parsimony table for the extant sequences
    """
    leaf_dict = {}
    for leaf in tree.iter_leaves():
        leaf_list = [(1, 0) if x == "-" else (0, 1) for x in leaf.sequence]
        leaf_dict[leaf.name] = leaf_list
    df = pd.DataFrame.from_dict(leaf_dict, orient="index")
    return df


def make_internal_parsimony(tree, extant_parsimony):
    """
    Make a parsimony table for the internal nodes
    :param tree: Tree with annotated sequences
    :param extant_parsimony: Parsimony table for the extant sequences
    :return: A parsimony table for the internal nodes
    """
    internal_dict = {}
    leaves = [x.name for x in tree.get_leaves()]

    for x in tree.traverse("postorder"):
        internal_list = []
        if not x.is_leaf():
            for column in extant_parsimony:
                child1_df = (
                    extant_parsimony
                    if x.children[0].name in leaves
                    else internal_parsimony
                )
                child2_df = (
                    extant_parsimony
                    if x.children[1].name in leaves
                    else internal_parsimony
                )

                child1 = child1_df.loc[x.children[0].name, column]
                child2 = child2_df.loc[x.children[1].name, column]

                internal_list.append(
                    (
                        min(child1[0], child1[1] + 1) + min(child2[0], child2[1] + 1),
                        min(child1[1], child1[0] + 1) + min(child2[1], child2[0] + 1),
                    )
                )

            internal_dict[x.name] = internal_list
        internal_parsimony = pd.DataFrame.from_dict(internal_dict, orient="index")

    return internal_parsimony


def get_phylo_contigs(parsimony_table):
    """
    Given a parsimony table return a mapping from node in the tree to which positions are present at that node
    :param tree:
    :param parsimony:
    :return:
    """

    phylo_contigs = defaultdict(list)
    filled_contigs = defaultdict(list)

    # Create a list of the columns where we have sequence content

    for row_idx, row in parsimony_table.iterrows():

        if row_idx != "N0":
            for col_idx, column in row.iteritems():
                if column[1] != 0:
                    phylo_contigs[row_idx].append(col_idx)

    # Fill in the contig dictionary - i.e. work out based on the sequence content positions where contiguous gaps /
    # sequence content exists

    print("phylo before")
    for k, v in phylo_contigs.items():
        print(k, v)

    for k, v in phylo_contigs.items():
        if k == "seq03":
            print("gotcha")
        filled_contigs[k] = fill_in_phylo_contigs(v, len(row))

    print("\n\nphylo filled")

    for k, v in filled_contigs.items():
        print(k, v)
    print()

    return filled_contigs


def get_path_contigs(phylo_contigs, path_marks, length, no_multiple_insertions=True):
    """Get a version of phylo contigs that is associated with a particular path and updated to reflect that there
    will be different points we can remove content from given the patterns we see in the ancestor"""

    path_contigs = defaultdict(lambda: defaultdict())

    print("here be hte path contigs")

    print(path_contigs)

    for path, vals in path_marks.items():
        for node, modes in vals.items():
            for mode, mark in sorted(modes.items()):
                if mode == "add":

                    path_contigs[path][node] = {"add": [mark]}

                if mode == "rmv":
                    if len(mark) > 1:

                        print("Mode is ", mode)
                        print("Node is ", node)
                        print("Path we're considering is ", path)
                        print("Pattern at extant is ", mark)

                        if node == "N1":
                            print("gotcha")

                        skip_pos = [x for x in range(length) if x not in path]

                        print("skip pos is ", skip_pos)
                        #
                        print("original phylo contigs is ", phylo_contigs)

                        contigs = [x for x in phylo_contigs[node]["add"]]

                        skipped_contigs = []

                        for contig in contigs:
                            for pos in contig:
                                if pos not in skip_pos:
                                    skipped_contigs.append(pos)

                        print("skipped contigs", skipped_contigs)

                        # for idx in modes['add']:
                        #     if idx not in skipped_contigs:
                        #         skipped_contigs.append(idx)

                        print("and now skipped contigs", skipped_contigs)

                        updated_phylo_contigs = fill_in_phylo_contigs(
                            skipped_contigs, length, skip=skip_pos
                        )

                        print("updated phylo contigs ", updated_phylo_contigs)

                        # If we aren't allowing for multiple insertions, remove any nodes that appear in the ancestor
                        #  from being added elsewhere
                        if no_multiple_insertions:
                            updated_phylo_contigs = remove_ancestral_insertions(
                                updated_phylo_contigs, path
                            )

                        if "add" in modes:

                            updated_phylo_contigs["add"] = [modes["add"]]

                        path_contigs[path][node] = updated_phylo_contigs

                        # for idx in modes['add']:
                        #     if idx not in updated_phylo_contigs[mode]:
                        #         updated_phylo_contigs[mode].append([idx])
                    else:
                        if node in path_contigs[path]:
                            path_contigs[path][node]["rmv"] = [mark]
                        else:
                            path_contigs[path][node] = {"rmv": [mark]}
    return path_contigs


def remove_ancestral_insertions(phylo_contigs, path):

    if path == (0, 2, 3, 5, 6, 7):
        print("gtocha")
    print("path is ", path)

    for mode, marks in phylo_contigs.items():
        if mode == "add":
            rmv = set()
            for idx_pos, idx in enumerate(marks):
                for pos_pos, pos in enumerate(idx):
                    if pos in path:
                        rmv.add(pos_pos)

                print("rmv was ", rmv)
                print(marks)
                for itm in sorted(rmv):
                    print(itm)
                    print(idx)
                    if len(idx) > itm:
                        idx.pop(itm)
                if len(idx) == 0:
                    marks.pop(idx_pos)

    return phylo_contigs


def get_highest_insertions(parsimony_table, tree):
    """
    Given a parsimony table, return the neccescary index for each column if inserting, to disallow multiple
    insertions in the same column
    :param parsimony_table:
    :return:
    """

    highest_insertions = {}

    for col in parsimony_table:
        insertion_candidates = []
        for node, score in parsimony_table[col].iteritems():

            if (
                node == "N0" and score[0] == 0
            ):  # If the root node has a zero cost for insertion
                highest_insertions[col] = node
                break
            else:
                if score[0] == 0:
                    insertion_candidates.append(node)
        if len(insertion_candidates) > 1:
            first = tree & insertion_candidates[0]
            others = [tree & x for x in insertion_candidates[1:]]
            common_ancestor = first.get_common_ancestor(others).name
            highest_insertions[col] = common_ancestor
        elif len(insertion_candidates) == 1:
            highest_insertions[col] = insertion_candidates[0]

    return highest_insertions


def make_total_parsimony(tree):
    """
    Make a combined extant and internal node parsimony table
    :param tree: Tree with annotated sequences
    :return: Combined parsimony table for extant sequences and internal nodes
    """

    extant_parsimony = make_extant_parsimony(tree)

    ancestor_parsimony = make_internal_parsimony(tree, extant_parsimony)

    frames = [extant_parsimony, ancestor_parsimony]

    total_parsimony = pd.concat(frames)

    return total_parsimony


def get_contig_indel_dict(parsimony_table):
    """
    Given a dataframe of parsimony values, return a dictionary mapping each row to the positions of any run of
    identical positions that are unbalanced between their parsimony scores."""

    indel_contigs = defaultdict(lambda: defaultdict(list))

    for row_idx, row in parsimony_table.iterrows():

        if row_idx != "N0":
            idx_list = []

            curr = (0, 0)

            for col_idx, column in row.iteritems():

                # If the cost here to add / remove a single column is zero, add it to the respective add / remove dictionary
                if column[0] == 0:
                    indel_contigs[row_idx]["add"].append([col_idx])
                elif column[1] == 0:
                    indel_contigs[row_idx]["rmv"].append([col_idx])

                # If it is the first column or it is identical to the previous column, add to the contiguous indels
                if col_idx == 0 or curr == column:
                    curr = column
                    idx_list.append(col_idx)

                # We've reached a column which doesn't have the same indel pattern, so add the previous columns to
                # indel_contigs
                if column != curr:

                    if len(idx_list) > 1:
                        if curr[0] == 0:
                            indel_contigs[row_idx]["add"].append(idx_list)

                        elif curr[1] == 0:
                            indel_contigs[row_idx]["rmv"].append(idx_list)

                    curr = column

                    idx_list = [col_idx]

                # If we've reached the end of the columns, add the index list to the contig indel dictionary
                if col_idx == len(parsimony_table.columns) - 1:

                    if len(idx_list) > 1:
                        if curr[0] == 0:
                            indel_contigs[row_idx]["add"].append(idx_list)

                        elif curr[1] == 0:
                            indel_contigs[row_idx]["rmv"].append(idx_list)

    return indel_contigs


def get_sibling(node):
    """
    Return the sibling of a node
    :param node: Node to find sibling of
    :return: Sibling if exists, else None
    """
    if node.up:
        children = node.up.get_children()
        sibling = [x for x in children if x.name != node.name]
        return sibling[0]
    else:
        return


def get_purged_contigs(contigs, tree, reduceContigs=True):
    """
    Given a contig dictionary and a tree, remove any optimal contigs that appear in higher up parts of the tree
    :param contigs: The contig dictionary
    :param tree: The tree annotated with sequences
    :param reduceContigs: Whether we should remove contigs that are a subset of other contigs
    :return: A possibly reduced contig dictionary that removes optimal contigs if they also appear at a higher level in
    the tree
    """
    purged_contigs = deepcopy(contigs)

    for node in tree.traverse("postorder"):

        # Can't collapse upwards into the root node
        if not node.is_root() and not node.up.is_root():

            sibling = get_sibling(tree.search_nodes(name=node.name)[0])

            # Purge for both add and remove
            for mode in purged_contigs[node.name].keys():

                for contig in contigs[node.name][mode]:

                    # If an identical contig occurs in a sibling, we can remove it from the current and sibling (as it
                    # could be completed at the parent node
                    if (
                        mode in purged_contigs[sibling.name]
                        and contig in purged_contigs[sibling.name][mode]
                    ):
                        purged_contigs[node.name][mode].remove(contig)
                        purged_contigs[sibling.name][mode].remove(contig)

                # If we want to remove contigs that are a subset of other contigs
                if reduceContigs:
                    purged_contigs[node.name][mode] = reduce_contigs(
                        purged_contigs[node.name][mode]
                    )

    return purged_contigs


def score_ancestor(paths, parsimony, contigs):
    path_scores = {}

    for path in paths:
        score = 0
        for col in parsimony:
            if col in path:
                score += parsimony.loc["N0"][col][0]
            else:
                score += parsimony.loc["N0"][col][1]
        path_scores[tuple(path)] = score
    return path_scores


def score_ancestor_collapse(paths, parsimony, contigs):
    path_scores = {}

    for path in paths:
        score = 0
        for col in parsimony:
            if col in path:
                score += parsimony.loc["N0"][col][0]
            else:
                score += parsimony.loc["N0"][col][1]
        path_scores[tuple(path)] = score
    return path_scores


def get_path_marks(paths, parsimony, contigs, aln, translate_path=True):

    path_scores = score_ancestor(paths, parsimony, contigs)
    path_marks = defaultdict(lambda: defaultdict(list))

    if translate_path:
        translation = alignment_permutations.get_translation(aln)

    for path, score in path_scores.items():
        if translate_path:
            translated_path = alignment_permutations.translate(path, translation)
        for col in parsimony:
            if (
                col in path
            ):  # If this is a column in our path we need to count the score for removing it
                score = parsimony.loc["N0"][col][0]
                if (
                    score != 0
                ):  # Score is non-zero, we need to account for removing content here
                    for k, v in contigs.items():
                        for event in v["rmv"]:
                            if len(event) > 0:
                                if col in event:

                                    if translate_path:
                                        path_marks[translated_path][k].append(
                                            "-%s-"
                                            % alignment_permutations.translate(
                                                [col], translation
                                            )
                                        )
                                    else:
                                        path_marks[path][k].append("-%s-" % str(col))

            else:  # This column isn't in our path so we need to count the score for adding it back in
                score = parsimony.loc["N0"][col][1]
                if score != 0:
                    for k, v in contigs.items():
                        for event in v["add"]:
                            if len(event) > 0:

                                if col in event:
                                    if translate_path:
                                        path_marks[translated_path][k].append(
                                            "+%s+"
                                            % alignment_permutations.translate(
                                                [col], translation
                                            )
                                        )
                                    else:
                                        path_marks[path][k].append("+%s+" % str(col))

    return path_marks


def get_path_marks2(
    paths,
    parsimony,
    contigs,
    aln,
    no_multple_insertions,
    highest_insertions,
    tree,
):

    path_scores = score_ancestor(paths, parsimony, contigs)
    path_marks = defaultdict(lambda: defaultdict(dict))

    for path, score in path_scores.items():

        for col in parsimony:
            if (
                col in path
            ):  # If this is a column in our path we need to count the score for removing it
                score = parsimony.loc["N0"][col][0]
                if (
                    score != 0
                ):  # Score is non-zero, we need to account for removing content here
                    for node, v in contigs.items():
                        if node not in path_marks[path]:
                            path_marks[path][node] = defaultdict(list)

                        # if 'rmv' not in path_marks[path][k]:
                        #
                        #     path_marks[path][k]["rmv"] = []

                        for event in v["rmv"]:
                            if len(event) > 0:

                                if (
                                    col in event
                                    and col not in path_marks[path][node]["rmv"]
                                ):

                                    path_marks[path][node]["rmv"].append(col)
                        if len(path_marks[path][node]) == 0:
                            path_marks[path].pop(node)

            else:  # This column isn't in our path so we need to count the score for adding it back in

                if (
                    no_multple_insertions
                ):  # We don't want to allow for inserting this position in multiple times

                    # print ("no multi")
                    # print (highest_insertions)

                    score = parsimony.loc["N0"][col][1]
                    # print ("score is ", score)
                    if score != 0:
                        for node, v in contigs.items():
                            if node not in path_marks[path]:

                                path_marks[path][node] = defaultdict(list)

                            # if 'add' not in path_marks[path][k]:
                            #     path_marks[path][k]["add"] = []
                            for event in v["add"]:

                                # If any of the nodes have a cost for adding this back in, record the node
                                if len(event) > 0:
                                    ins_event = highest_insertions[col]
                                    # print (ins_event)
                                    break
                                    # if col in event and col not in path_marks[path][node]['add']:
                                    #     path_marks[path][node]["add"].append(col)

                        if "add" not in path_marks[path][ins_event]:
                            path_marks[path][ins_event]["add"] = []

                        # Add the insertion event at the highest neccescary point to stop multiple insertions
                        path_marks[path][ins_event]["add"].append(col)

                        if len(path_marks[path][node]) == 0:
                            path_marks[path].pop(node)

                            # Now we need to check if any descendants need to have this insertion removed
                            curr = tree & ins_event

                            descendants = [x.name for x in curr.get_descendants()]

                            for desc in descendants:
                                for pos in contigs[desc]["rmv"]:
                                    if col in pos:
                                        if "rmv" not in path_marks[path][desc]:
                                            path_marks[path][desc]["rmv"] = []

                                        if col not in path_marks[path][desc]["rmv"]:
                                            print("adding")
                                            path_marks[path][desc]["rmv"].append(col)
                            print()

                else:  # We will allow to insert this position in multiple times

                    score = parsimony.loc["N0"][col][1]
                    if score != 0:
                        for node, v in contigs.items():
                            if node not in path_marks[path]:

                                path_marks[path][node] = defaultdict(list)

                            # if 'add' not in path_marks[path][k]:
                            #     path_marks[path][k]["add"] = []
                            for event in v["add"]:
                                # If any of the nodes have a cost for adding this back in, record the node
                                if len(event) > 0:
                                    if (
                                        col in event
                                        and col not in path_marks[path][node]["add"]
                                    ):
                                        path_marks[path][node]["add"].append(col)
                            if len(path_marks[path][node]) == 0:
                                path_marks[path].pop(node)

    # print ('path marks here is ')
    # for k,v in path_marks.items():
    #     for x, y in v.items():
    #         print (k, x, y)
    return path_marks


def concat_lists(marks, contigs, phylo_contigs):
    for contig in contigs:
        for idx, col in enumerate(contig):
            if col not in marks:
                contig.pop(idx)
    return contigs[0]


def concatenate_marks(path_marks, contigs, phylo_contigs):
    # Check to see if we can concatenate any of the path marks
    concatenated_marks = deepcopy(path_marks)

    for path, marks in path_marks.items():

        for node, mark in marks.items():
            if path == (0, 4):
                break

                for mode in mark:
                    if len(mark[mode]) > 1:
                        concatenated = concat_lists(
                            mark[mode], contigs[node][mode], [0, 1]
                        )
                        concatenated_marks[path][node][mode] = concatenated
                    else:
                        pass
                        # # print ("mark", mark)
                        # # print ("contigs[node]",contigs[node])
                        # # print (contigs[node].values())
                        # for val in contigs[node][mode]:
                        #     # print ('val is ', val)
                        #     if mark == val:
                        #         pass
                        #         # print ("HERE")
                        #     # if translate_path:
                        #     #     concatenated_marks[translated_path][k].append(
                        #     #         "+%s+" % alignment_permutations.translate([col], translation))
                        #     else:
                        #         # print ("got here")
                        #         concatenated_marks[path][node] = []
                        #         concatenated_marks[path][node].append(",".join(str(x) for x in val))
                        #         break
                        # concatenated_marks[path][k].append("SOMETHING")

    return concatenated_marks


def get_marks_to_ancestors(path_marks):
    """
    Return a dictionary mapping, for each mode, a given node to the positions that can be dealt with there
    :return: The dictionary
    """

    marks_to_ancestors = defaultdict(lambda: defaultdict(lambda: (defaultdict(list))))

    for path, marks in path_marks.items():
        for node, mark in marks.items():
            for mode in mark:
                marks_to_ancestors[node][mode][path].append(mark[mode])

                # marks_to_ancestors[node][mode].add({path: mark[mode]})

    return marks_to_ancestors


def shift_marks(path_marks, contigs, tree):

    shifted_marks = deepcopy(path_marks)
    marks_to_ancestors = get_marks_to_ancestors(path_marks)

    for node in tree.traverse("preorder"):
        if not node.is_root() and not node.is_leaf():

            children = node.children

            for mode in marks_to_ancestors[node.name]:

                if (
                    mode in marks_to_ancestors[children[0].name]
                    and mode in marks_to_ancestors[children[1].name]
                ):
                    for path in marks_to_ancestors[node.name][mode]:
                        if (
                            path in marks_to_ancestors[children[0].name][mode]
                            and path in marks_to_ancestors[children[1].name][mode]
                        ):
                            for idx in marks_to_ancestors[node.name][mode][path]:
                                for pos in idx:
                                    present = []
                                    for idx_child1 in marks_to_ancestors[
                                        children[0].name
                                    ][mode][path]:
                                        for pos_child1 in idx_child1:
                                            if pos == pos_child1:
                                                present.append(children[0].name)
                                    for idx_child2 in marks_to_ancestors[
                                        children[1].name
                                    ][mode][path]:
                                        for pos_child2 in idx_child2:
                                            if pos == pos_child2:
                                                present.append(children[1].name)
                                    # print ('here is present ', present)
                                    for child_name in present:
                                        if pos in shifted_marks[path][child_name][mode]:

                                            shifted_marks[path][child_name][
                                                mode
                                            ].remove(pos)
                                            # for subchildren in tree&child_name:
                                            #     print ('here is a subchild')
                                            #     print (subchildren.name)

    return shifted_marks


def move_marks(path_marks, path_contigs, tree, length, push_down=False):

    shifted_marks = deepcopy(path_marks)
    marks_to_ancestors = get_marks_to_ancestors(path_marks)
    pushed_down = defaultdict(list)

    for node in tree.traverse("postorder"):
        if not node.is_root() and not node.is_leaf():

            children = node.children

            for mode in marks_to_ancestors[node.name]:

                if node.name == "N1":
                    print("gotcha")

                if (
                    mode in marks_to_ancestors[children[0].name]
                    and mode in marks_to_ancestors[children[1].name]
                ):
                    for path in marks_to_ancestors[node.name][mode]:
                        if (
                            path in marks_to_ancestors[children[0].name][mode]
                            and path in marks_to_ancestors[children[1].name][mode]
                        ):
                            # print ('no sno ', path_contigs[path])
                            # print (node.name)
                            # print ('here bo ', path_contigs[path][node.name])
                            if mode in path_contigs[path][node.name]:
                                for idx in path_contigs[path][node.name][mode]:

                                    # for idx in path_marks[path][node.name][mode]:
                                    # for idx in marks_to_ancestors[node.name][mode][path]:

                                    for pos in idx:
                                        present = []
                                        found_idxes = []
                                        for idx_child1 in marks_to_ancestors[
                                            children[0].name
                                        ][mode][path]:
                                            for pos_child1 in idx_child1:
                                                if pos == pos_child1:
                                                    present.append(children[0].name)
                                                    found_idxes.append(idx_child1)
                                                    break
                                        for idx_child2 in marks_to_ancestors[
                                            children[1].name
                                        ][mode][path]:
                                            for pos_child2 in idx_child2:
                                                if pos == pos_child2:
                                                    present.append(children[1].name)
                                                    found_idxes.append(idx_child2)
                                                    break

                                        # print ('here is present ', present)
                                        # if len(idx) == 1 or (found_idxes[0] == found_idxes[1]):

                                        if len(found_idxes) > 1:

                                            if found_idxes[0] == found_idxes[1]:

                                                print("triple doo")
                                                for child_name in present:
                                                    if (
                                                        pos
                                                        in shifted_marks[path][
                                                            child_name
                                                        ][mode]
                                                    ):
                                                        shifted_marks[path][child_name][
                                                            mode
                                                        ].remove(pos)
                                                        # for subchildren in tree&child_name:
                                                        #     print ('here is a subchild')
                                                        #     print (subchildren.name)
                                            elif found_idxes[0] != found_idxes[1]:
                                                print("woo poo")
                                                if (
                                                    pos
                                                    in shifted_marks[path][node.name][
                                                        mode
                                                    ]
                                                ):
                                                    shifted_marks[path][node.name][
                                                        mode
                                                    ].remove(pos)
                                                pushed_down[pos].append(present[0])
                                                pushed_down[pos].append(present[1])
                                                print(
                                                    "pushed down is ",
                                                    pushed_down,
                                                )

                                                # If we are removing something that was pushed down earlier,
                                                # revert the push down operation
                                                if push_down:
                                                    if pos in pushed_down:
                                                        for seq in pushed_down[pos]:
                                                            print(
                                                                node.get_descendants()
                                                            )
                                                            if seq in [
                                                                x.name
                                                                for x in node.get_descendants()
                                                            ]:
                                                                print(
                                                                    shifted_marks[path][
                                                                        seq
                                                                    ]
                                                                )
                                                                if (
                                                                    pos
                                                                    not in shifted_marks[
                                                                        path
                                                                    ][
                                                                        seq
                                                                    ][
                                                                        mode
                                                                    ]
                                                                ):
                                                                    shifted_marks[path][
                                                                        seq
                                                                    ][mode].append(pos)

                                                if not push_down:
                                                    if pos in pushed_down:
                                                        shifted_marks[path][node.name][
                                                            mode
                                                        ].append(pos)
                                                        for seq in pushed_down[pos]:
                                                            print(
                                                                node.get_descendants()
                                                            )
                                                            if seq in [
                                                                x.name
                                                                for x in node.get_descendants()
                                                            ]:
                                                                print(
                                                                    shifted_marks[path][
                                                                        seq
                                                                    ]
                                                                )
                                                                if (
                                                                    pos
                                                                    in shifted_marks[
                                                                        path
                                                                    ][seq][mode]
                                                                ):
                                                                    shifted_marks[path][
                                                                        seq
                                                                    ][mode].remove(pos)

                                                print(idx)
                                                print(found_idxes)
                                                print("hmm")
    return shifted_marks


def separate_contig(marks, phylo_contig):
    separated_contigs = []
    new_contig = []

    for contigs in phylo_contig:
        for mark in marks:
            if mark in contigs:
                new_contig.append(mark)
            elif len(new_contig) > 0:
                separated_contigs.append(new_contig)
                new_contig = []
    if len(new_contig) > 0 and new_contig not in separated_contigs:
        separated_contigs.append(new_contig)
    return separated_contigs


def get_sequence_content_from_ancestors(path_marks, node, tree):
    # print ('inthis one')
    # print ('here is path marks')
    # print (path_marks)
    #
    # print (node)
    # print (tree&node)
    sequence_content = []
    tree_node = tree & node
    for ancestor in tree_node.get_ancestors():
        # print ("ANCESTOR IS ", ancestor.name)
        if ancestor.name in path_marks:
            sequence_content.append(x for x in path_marks[ancestor.name]["add"])
            # print ("her eis ", path_marks[ancestor.name])
    # print (tree_node.get_ancestors())

    # print ("sequence content is ", sequence_content)

    return sequence_content


def separate_marks(path_marks, phylo_contigs, tree, length):

    separated_marks = deepcopy(path_marks)

    for path, vals in path_marks.items():
        for node, modes in vals.items():
            for mode, mark in modes.items():
                if mode == "rmv":
                    if len(mark) > 1:
                        # print ('marks is ', mark)
                        # if path == (0,4) and node == "N4":
                        #     print ('gotcha')
                        # if path == (0,1,3):
                        #     print ('got this one')

                        #
                        print("Mode is ", mode)
                        print("Node is ", node)
                        print("Path we're considering is ", path)
                        print("Pattern at extant is ", mark)

                        # Get the updated sequence content at this node (relative to what has been placed in its ancestors)
                        # seq_content_from_ancestor = get_sequence_content_from_ancestors(path_marks[path],  node, tree)

                        # Add in another that is identical between the root ancestor and the extant

                        # seq_content_from_ancestor += [x for x in path]
                        # print ("Seq content from ancestor is ", seq_content_from_ancestor)
                        #
                        # print ('length here is ', length)

                        skip_pos = [x for x in range(length) if x not in path]

                        print("skip pos is ", skip_pos)
                        #
                        print("original phylo contigs is ", phylo_contigs)

                        contigs = [x for x in phylo_contigs[path][node]["add"]]

                        # print ('now here be the skipped')

                        # print (contigs)

                        skipped_contigs = []

                        for contig in contigs:
                            for pos in contig:
                                if pos not in skip_pos:
                                    skipped_contigs.append(pos)

                        print("skipped contigs", skipped_contigs)

                        updated_phylo_contigs = fill_in_phylo_contigs(
                            skipped_contigs, length, skip=skip_pos
                        )

                        print("updated phylo contigs ", updated_phylo_contigs)

                        # print (updated_phylo_contigs)

                        # Get the cost for placing
                        separated_contigs = separate_contig(
                            mark, updated_phylo_contigs[mode]
                        )

                        print("separated contigs ", separated_contigs)

                        # separated_contigs = separate_contig(mark, phylo_contigs[node][mode])
                        separated_marks[path][node][mode] = separated_contigs

                        print("separated marks ", separated_marks)

    return separated_marks


def get_event_count(path):
    """
    Count the number of indel events implied by a particular ancestral path
    :param path: The ancestral path
    :return: The number of indel events
    """

    events = [x for x in path.values()]
    count = 0
    for events in path.values():
        event = "".join(x for x in events)
        count += (event.count("+") + event.count("-")) / 2

    return int(count)


def translate_paths2(path_marks, aln, translation=None):

    translated_paths = defaultdict(lambda: defaultdict(list))

    if not translation:
        translation = alignment_permutations.get_translation(aln)

    for path, nodes in path_marks.items():

        translated_path = alignment_permutations.translate("path", path, translation)
        for node, pos in nodes.items():
            # if path == (0,1,2,3,4):
            #     print ('gotcha')
            for mode, val in pos.items():
                if len(val) > 0:

                    if len(mode) > 0:
                        translated_node = alignment_permutations.translate(
                            mode, val, translation
                        )
                        translated_paths[translated_path][node].append(translated_node)

    return translated_paths


def get_marked_style(title, path):
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_leaf_name = False
    ts.branch_vertical_margin = 10  # 10 pixels between adjacent branches
    event_count = get_event_count(path)

    ts.title.add_face(
        TextFace(
            "Minimal indel events needed for %s ancestor is %s " % (title, event_count),
            fsize=20,
        ),
        column=0,
    )

    def my_layout(node):
        if node.name == "N0":  # Format the root node with the ancestor sequence
            F = TextFace(title, tight_text=True, fgcolor="black")
            # F.rotation = 330
            add_face_to_node(F, node, column=0, position="branch-top")

        if node.name in path.keys():

            annot = "".join(path.get(node.name))
            F = TextFace(annot, tight_text=True, fgcolor="red")
            F.rotation = 330
            add_face_to_node(F, node, column=0, position="branch-bottom")
        if node.is_leaf():
            seq_face = SeqMotifFace(node.sequence, seqtype="aa", seq_format="seq")
            add_face_to_node(seq_face, node, column=0, position="aligned")

    ts.layout_fn = my_layout

    return ts


def reduce_contigs(contig_list):
    """
    Method to remove any contig position in a contig list if it is a subset of another contig
    """
    reduced_list = contig_list[:]
    for contig in contig_list:
        for other_contig in reduced_list:
            if set(contig).issubset(set(other_contig)) and contig != other_contig:
                reduced_list.remove(contig)
                break
    return reduced_list


def is_next(pos, check_pos):
    """
    Check if two positions are contiguous. Currently defined as being +1
    """
    if check_pos == pos + 1:
        return True
    else:
        return False


def split_contigs(query_contig):
    """
    Given a list of contiguous positions, reduce the list to a new list of lists of the actual contigs
    """
    new_contigs = []
    start_pos = 0
    for pos, x in enumerate(zip(query_contig, query_contig[1:])):
        if not is_next(x[0], x[1]):
            new_contigs.append(query_contig[start_pos : pos + 1])
            start_pos = pos + 1

    # Add the final segment
    new_contigs.append(query_contig[start_pos:])

    return new_contigs


def produce_split_contig_list(contig_list):
    split_contig_list = []
    for contig in contig_list:
        for new_contig in split_contigs(contig):
            if len(new_contig) > 0:
                split_contig_list.append(new_contig)
    return split_contig_list


def remove_positions(query, target):
    """
    Remove all of the individual positions within a list of lists that don't appear in a list of queries"""

    for sublist in target:
        for pos in sublist:
            if pos not in query:
                sublist.remove(pos)

    return target


class IndelPlacement:
    def __init__(
        self,
        tree_name,
        aln_name,
        get_all_ancestors=False,
        no_multiple_insertions=True,
    ):
        self.tree_name = tree_name
        self.aln_name = aln_name
        self.get_all_ancestors = get_all_ancestors

        self.tree = PhyloTree(
            tree_name, alignment=aln_name, format=1, alg_format="fasta"
        )
        self.aln = AlignIO.read(aln_name, "fasta")

        label_internal_nodes(self.tree)

        self.parsimony = make_total_parsimony(self.tree)

        self.parsimony.sort_index(inplace=True)

        self.phylo_contigs = get_phylo_contigs(self.parsimony)

        self.contigs = get_contig_indel_dict(self.parsimony)

        self.highest_insertions = get_highest_insertions(self.parsimony, self.tree)

        # Change this back
        self.purged_contigs = get_contig_indel_dict(self.parsimony)

        # self.purged_contigs = get_purged_contigs(self.contigs, self.tree)

        if self.get_all_ancestors:
            self.paths = alignment_permutations.get_all_paths(self.aln)

        else:
            self.paths = alignment_permutations.get_ancestral_paths(self.aln)

        insertion_aware_ancestral_paths = []
        if no_multiple_insertions:
            ancestor_positions = [
                k for k, v in self.highest_insertions.items() if v == "N0"
            ]
            # print(ancestor_positions)
            for path in self.paths:
                if set(ancestor_positions).issubset(set(path)):
                    # print(path)
                    # print("was subset")
                    insertion_aware_ancestral_paths.append(path)

            self.paths = insertion_aware_ancestral_paths

        # print ('here da paths')
        # print (self.paths)

        # self.paths = [[0, 1, 4, 5]]

        self.scores = score_ancestor(self.paths, self.parsimony, self.purged_contigs)

        self.marks = get_path_marks2(
            self.paths,
            self.parsimony,
            self.purged_contigs,
            self.aln,
            no_multiple_insertions,
            self.highest_insertions,
            self.tree,
        )

        self.path_contigs = get_path_contigs(
            self.phylo_contigs,
            self.marks,
            len(self.parsimony.columns),
            no_multiple_insertions,
        )

        print("\n\n\n")
        print("phylo contigs")

        for k, v in self.phylo_contigs.items():
            for a, b in v.items():
                print(k, a, b)

        print("\n")

        print("path contigs")

        for k, v in self.path_contigs.items():
            for a, b in v.items():
                print(k, a, b)
        print("\n\n\n")

        # NOTE: NO LONGER RUNNING CONCATENATE

        # self.concatenated_marks = concatenate_marks(self.marks, self.purged_contigs, self.phylo_contigs)

        # self.shifted_marks = shift_marks(self.marks, self.purged_contigs, self.tree)

        self.moved_marks = move_marks(
            self.marks,
            self.path_contigs,
            self.tree,
            len(self.parsimony.columns),
        )

        print(" marks")
        for k, v in self.marks.items():
            for a, b in v.items():
                print(k, a, b)

        # print ('concatenated marks')
        # for k, v in self.concatenated_marks.items():
        #     if k == (0,3,4):
        #         for a, b in v.items():
        #             print (k,a, b)

        # print ('shifted marks')
        # for k, v in self.shifted_marks.items():
        #     if k == (0,3,4):
        #         for a, b in v.items():
        #             print (k,a, b)

        print("moved marks")
        for k, v in self.moved_marks.items():
            for a, b in v.items():
                print(k, a, b)

        self.separated_marks = separate_marks(
            self.moved_marks,
            self.path_contigs,
            self.tree,
            len(self.parsimony.columns),
        )

        print("separated marks")

        for k, v in self.separated_marks.items():
            for a, b in v.items():
                print(k, a, b)

        self.translated_marks = translate_paths2(self.separated_marks, self.aln)

        # print ('translated marks')
        # for k, v in self.translated_marks.items():
        #     print (k, v)
