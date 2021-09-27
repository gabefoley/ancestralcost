import itertools
import ancestralcost.graph as graph


def paths(graph, v):
    """Generate the maximal cycle-free paths in graph starting at v.
    graph must be a mapping from vertices to collections of
    neighbouring vertices.

    # >>> g = {1: [2, 3], 2: [3, 4], 3: [1], 4: []}
    # >>> sorted(paths(g, 1))
    # [[1, 2, 3], [1, 2, 4], [1, 3]]
    # >>> sorted(paths(g, 3))
    [[3, 1, 2, 4]]

    """
    path = [v]  # path traversed so far
    seen = {v}  # set of vertices in path

    def search():
        dead_end = True
        for neighbour in graph[path[-1]]:
            if neighbour not in seen:
                dead_end = False
                seen.add(neighbour)
                path.append(neighbour)
                yield from search()
                path.pop()
                seen.remove(neighbour)
        if dead_end:
            yield list(path)

    yield from search()


def get_unique_ancestors_from_aln(aln):
    """
    Given the pattern seen in an alignment file, generate a list of possible sequences, with a * added as a stop
    :return: List of sequences with * in the final position of each sequence
    """
    unique_seqs = []
    for seq in aln:
        modified_seq = str(seq.seq) + "*"
        if modified_seq not in unique_seqs:
            unique_seqs.append(modified_seq)

    return unique_seqs


def translate(mode, index, translation):

    if isinstance(index, tuple):
        return "".join([translation[int(x)] for x in index])
    elif type(index) == list:
        translated = []

        symbol = "+" if mode == "add" else "-"

        if type(index[0]) == list:
            translation_string = ""
            for i in index:
                translation_string += (
                    symbol + "".join([translation[int(x)] for x in i]) + symbol + " "
                )
                # translated.append(translation_string)

        else:
            translation_string = (
                symbol + "".join([translation[int(x)] for x in index]) + symbol
            )

        translated.append(translation_string.strip())

        return translated[0]

        # indel_type = "+" if "+" in index[0] else "-"
        #
        # print("INDEX IS")
        # print(index[0])
        # print("INDEL TYPE IS ")
        # print(indel_type)
        # print ("Stripped version is ")
        # m = re.match(r"[+](.*)\.[+]", index[0])
        # print ("Regex match is")
        # print (m.group(1))
        # replaced = [x.replace(indel_type, "") for x in index][0]
        # print (replaced)

        # return indel_type + "".join([translation[int(x.strip(indel_type))] for x in replaced]) + indel_type


def get_translation(aln):

    translation = {}
    translation[-1] = ""
    translation[len(aln[0])] = "*"

    for seq in range(len(aln)):
        for pos, sym in enumerate(aln[seq]):
            if pos not in translation and sym != "-":
                translation[pos] = sym

    return translation


def get_traversals(adj_dict, start_pos=-1):
    """
    Given an adjacency dictionary, get a list of all the possible path
    :param adj_dict: Adjacency dictionary
    :param start_pos: Index to start from, default is -1
    :return: List of all the traversals
    """
    return sorted(paths(adj_dict, start_pos))


def get_ancestral_paths(aln, single=False):
    """
    Get a list of all the valid ancestral paths as implied by the alignment
    :param aln: The alignment
    :return: A list of all valid ancestral paths
    """

    unique_seqs = get_unique_ancestors_from_aln(aln)

    if single:
        adj_dict = graph.get_adjaceny_dict([unique_seqs[0]])

    else:
        adj_dict = graph.get_adjaceny_dict(unique_seqs)

    traversals = get_traversals(adj_dict)

    paths = []

    for path in traversals:
        paths.append(path[1:-1])

    return paths


def get_all_paths(aln):
    """
    From an alignment get all possible combinations of ancestors (assumes linear ordering, but sequences can start or
    stop at any point)
    :param aln: The alignment
    :return: List of all possible combinations of ancestors
    """
    seq_list = []

    numbers = [x for x in range(len(aln[0]))]

    for i in range(len(numbers) + 1):
        for item in itertools.combinations(numbers, i):

            path = [int(x) for x in item]
            if path not in seq_list:
                seq_list.append(path)

    return seq_list
