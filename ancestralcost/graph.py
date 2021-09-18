# Python program to print topological sorting of a DAG
from collections import defaultdict


# Class to represent a graph
class Graph:
    def __init__(self, vertices):
        self.graph = defaultdict(list)  # dictionary containing adjacency List
        self.V = vertices  # No. of vertices

    # function to add an edge to graph
    def addEdge(self, u, v):
        self.graph[u].append(v)

        # A recursive function used by topologicalSort

    def topologicalSortUtil(self, v, visited, stack):

        # Mark the current node as visited.
        visited[v] = True

        # Recur for all the vertices adjacent to this vertex
        for i in self.graph[v]:
            if visited[i] is False:
                self.topologicalSortUtil(i, visited, stack)

                # Push current vertex to stack which stores result
        stack.insert(0, v)

        # The function to do Topological Sort. It uses recursive

    # topologicalSortUtil()
    def topologicalSort(self):
        # Mark all the vertices as not visited
        visited = [False] * self.V
        stack = []

        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
        for i in range(self.V):
            if visited[i] is False:
                self.topologicalSortUtil(i, visited, stack)

                # Print contents of the stack
        # print (stack)


def get_adjaceny_dict(seqs):

    adj_dict = defaultdict(set)

    for seq in seqs:
        prev_node = -1
        for pos, sym in enumerate(seq):
            curr_node = pos
            if sym != "-":
                adj_dict[prev_node].add(curr_node)
                prev_node = pos

    return adj_dict


def make_graph(seqs):
    adj_dict = get_adjaceny_dict(seqs)

    graph = Graph(len(seqs[0]))

    for node, next_nodes in adj_dict.items():
        for next_node in next_nodes:
            graph.addEdge(node, next_node)
    return graph
