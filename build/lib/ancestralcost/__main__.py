""" Check for each position in a given ancestor that the presence of ancestral content implied to be there by a given alignment and tree is not substantially less parsimonious then the alternative of not having ancestral content there."""

import argparse
import sys

from ete3 import PhyloTree
import ancestralcost.indel_placement as indel_placement
import ancestralcost.tree_code as tree_code

from Bio import AlignIO

from ancestralcost import __version__

def run(args):
    print(f"\U0001F4B0" + " Running Ancestral Cost v" + __version__ + " \U0001F4B0" + "\n")

    node_of_interest = "#N0" if not args.node else args.node


    aln = AlignIO.read(args.aln, "fasta")

    tree = tree_code.load_tree_with_alignment(args.tree, args.aln)

    if args.fasta_output:
        print("Getting the positions required to be there for all ancestors")
        with open(args.fasta_output, "w+") as fasta_output:
            for node in tree.traverse():
                if not node.is_leaf():
                    print(node.name)
                    positions = tree_code.get_positions_with_content(tree, aln, node.name)
                    print(positions)
                    print(len(aln[0]))

                    ancestor_gaps = [
                        "A" if x in positions else "-" for x in range(len(aln[0]))
                    ]

                    print(ancestor_gaps)
                    fasta_output.write(f'>{node.name}\n{"".join(ancestor_gaps)}\n')

        with open(args.tree_output, "w+") as tree_output:
            print(tree.write(format=7))
            tree_output.write(tree.write(format=1))

    else:

        print("Getting the positions that are required to be there\n")
        positions = tree_code.get_positions_with_content(tree, aln, node_of_interest)

        print(
            "These positions are required to be there - ",
            [x + 1 for x in positions],
        )

        # If we just want the positions, print them here and don't continue the rest of the progam
        if not args.positions:

            print("\nLabelling internal nodes\n")
            previous_labels = indel_placement.label_internal_nodes(tree)

            node_of_interest = (
                node_of_interest if not args.node else previous_labels[args.node]
            )

            print("Getting the parsimony scores\n")
            parsimony = tree_code.get_parsimony_table(tree)

            print("Getting the parsimony scores for required positions\n")
            parsimony_scores = tree_code.get_parsimony_scores(
                parsimony, positions, node_of_interest
            )

            print("Retrieving values for " + node_of_interest + "\n")

            if not parsimony_scores:
                print("No positions had unbalanced parsimony scores")

            else:

                for (k, v) in parsimony_scores.items():
                    print(
                        f"Position {k + 1} was required to be in the ancestor,"
                        " but the parsimony cost for it being there was"
                        f" {v[0]} compared to the parsimony cost for it not"
                        f" being there which was {v[1]}"
                    )
                    print()



def ac_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--aln", help="Path to alignment", required=True)
    parser.add_argument("-t", "--tree", help="Path to phylogenetic tree", required=True)
    parser.add_argument(
        "-n", "--node", help="Node to return cost for (default is root)"
    )
    parser.add_argument(
        "-p",
        "--positions",
        help="Just return the positions required to be there",
        action="store_true",
    )
    parser.add_argument(
        "-f", "--fasta_output", help="Return all ancestors as a FASTA file"
    )
    parser.add_argument("-to", "--tree_output", help="Write out the ancestor tree")

    return parser


def main(args=None):
    parser = ac_parser()
    if args:
        sys.argv = args
    if len(sys.argv) > 1 and sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f'ancestralcost v{__version__}')
        sys.exit(0)
    else:
        print(f'ancestralcost v{__version__}')
        args = parser.parse_args(args)
        run(args)
    sys.exit(0)


if __name__ == "__main__":
    main()