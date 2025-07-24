from io import StringIO
from Bio import Phylo


class NeighborJoining:
    def __init__(self, dist_matrix):
        self.dist_matrix = dist_matrix

    # // distance calculation stuff for NJ
    # returns a list of total distances for forming Distance Matrix Prime
    def _total_row_distances(self):
        n = len(self.dist_matrix)
        total_distances = []
        for i in range(n):
            total_sum = 0
            for j in range(n):
                total_sum += self.dist_matrix[i][j]
            total_distances.append(total_sum)
        return total_distances

    # adjustedDistanceFollows a different calculation instead
    def _adjusted_distance(self, divergences):
        adj_matrix = []
        mat_len = len(self.dist_matrix)
        for i in range(mat_len):
            node_row = []
            adj_matrix.append(node_row)
            for j in range(mat_len):
                if j == i:
                    node_row.append(0)
                else:
                    dIJ = self.dist_matrix[i][j]
                    dIJ_prime = ((mat_len - 2) * dIJ) - (
                        divergences[i] + divergences[j]
                    )
                    node_row.append(dIJ_prime)
        return adj_matrix

    def _pair_distance(self, nodeI, nodeJ):
        # return new calculated distances
        stored_values = []
        mat_len = len(self.dist_matrix)
        for k in range(mat_len):
            if k != nodeI and k != nodeJ:
                # dMI/dMJ
                dM = (
                    self.dist_matrix[nodeI][k]
                    + self.dist_matrix[nodeJ][k]
                    - self.dist_matrix[nodeI][nodeJ]
                ) / 2
                stored_values.append(dM)

        return stored_values

    # limb length is calculated slightly difference by taking the delta between
    # nodes A and B into consideration instead of divergences
    # calculate limb lengths for each leaf that is joined
    # return a tuple containing two values for each distance
    def _limb_length(self, nodeA, nodeB, divergences):
        n = len(self.dist_matrix)
        dAB = self.dist_matrix[nodeA][nodeB]
        divergenceA = divergences[nodeA]
        divergenceB = divergences[nodeB]
        deltaAB = (divergenceA - divergenceB) / (n - 2)
        # limb lengths
        dAZ = (dAB + deltaAB) / 2
        dBZ = (dAB - deltaAB) / 2

        return dAZ, dBZ

        # // Cluster NJ Stuff

    def _to_newick(self, tree):
        def recurse(node):
            if isinstance(tree[node], dict):
                children = []
                for child, dist in tree[node].items():
                    if child in tree:
                        children.append(f"{recurse(child)}:{dist}")
                    else:
                        children.append(f"{child}:{dist}")
                return f"({','.join(children)})"
            return node

        # Get the topmost node (root)
        root = list(tree.keys())[-1]
        return recurse(root) + ";"

    def _cluster_NJ(self, tree, nodes):
        mat_len = len(self.dist_matrix)
        if mat_len == 2:
            return self.dist_matrix
        divergences = self._total_row_distances()
        adj_distance_matrix = self._adjusted_distance(divergences)
        min_val = float("inf")
        min_i, min_j = 0, 0
        # looking for neighboring nodes,  based on minimum value, node i and node j
        # will be neighbors
        for i in range(mat_len):
            for j in range(mat_len):
                val = adj_distance_matrix[i][j]
                if val < min_val:
                    min_val = val
                    # store node indices
                    min_i = i
                    min_j = j
        # merge node I and node J together
        # grab distances of other nodes
        new_limbs = self._limb_length(min_i, min_j, divergences)
        # these limbs are for tree construction
        new_limb_MI, new_limb_MJ = new_limbs[0], new_limbs[1]
        # styling choice for tree
        tree[f"({nodes[min_j]}<>{nodes[min_i]})"] = {
            nodes[min_i]: new_limb_MI,
            nodes[min_j]: new_limb_MJ,
        }
        """
        NEED TO FIX THIS LINE BEFORE MERGING
        currently feng_doolittle only works with 2 or 3 sequences
        single sequences or >3 sequences have index error
        """
        nodes.insert(nodes.index(nodes[min_i]), f"({nodes[min_j]}<>{nodes[min_i]})")
        nodes.remove(str(nodes[min_i + 1]))
        nodes.remove(str(nodes[min_j]))

        # calculate new distances for new node to remaining nodes
        # construct a new distance matrix
        # store distance matrice values for new node
        new_node_distances = self._pair_distance(min_i, min_j)

        # filler node that will be not used, only there to skip iteration
        new_distance_matrix = [
            [0 for _ in range(mat_len - 1)] for _ in range(mat_len - 1)
        ]
        saved_matrices_values = []

        # k is the row_idx and m is the col_idx, change if desired
        # grabbing values that have not been merged
        for i in range(len(new_node_distances)):
            saved_matrices_values.append(new_node_distances[i])
        for k in range(mat_len):
            if k == min_i or k == min_j:
                continue
            for m in range(mat_len - 1):
                if m == min_i or m == min_j:
                    continue
                if self.dist_matrix[k][m] != 0:
                    saved_matrices_values.append(self.dist_matrix[k][m])
        if mat_len - 1 > 2:
            # Fill in newDistanceMatrix based on conditions
            count = len(saved_matrices_values) - 1
            for i in range(mat_len):
                for j in range(i + 1, mat_len):
                    if j - 1 == i:
                        new_distance_matrix[i][j - 1] = 0
                        new_distance_matrix[j - 1][i] = 0
                    else:
                        new_distance_matrix[i][j - 1] = saved_matrices_values[j - count]
                        new_distance_matrix[j - 1][i] = saved_matrices_values[j - count]
                count -= 1
                # now solve i for when j does not equal zero
        elif mat_len - 1 == 2:
            new_distance_matrix[0][1] = saved_matrices_values[0]
            new_distance_matrix[1][0] = saved_matrices_values[0]
        # print(new_distance_matrix)
        return new_distance_matrix, tree, nodes

    # place holder NJ algorithm formatting to Newick
    def generate_newick(self):
        # distance matrix and n x n dimensions, respectfully
        dist_matrix = self.dist_matrix
        tree = {}
        nj_nodes = [str(i) for i in range(len(dist_matrix))]
        while len(dist_matrix) != 2:
            dist_matrix, tree, nj_nodes = self._cluster_NJ(tree, nj_nodes)
            # merge remaining nodes in 2x2
        # perform merge on final 2 nodes
        node_a, node_b = nj_nodes[0], nj_nodes[1]
        dist = dist_matrix[0][1]
        limb_length = dist / 2

        # Make sure internal node is formatted correctly
        last_key = f"{node_a}<>{node_b}"
        tree[last_key] = {node_a: limb_length, node_b: limb_length}
        # clean up tree so nodes reflect <> notation
        final_tree = {}
        for key in tree:
            if key[0] == "[":
                final_tree[key[1 : len(key) - 1]] = tree[key]
            else:
                final_tree[key] = tree[key]
        return self._to_newick(final_tree)


# takes a matrix from a clustering algorithm and outputs a newick tree, can also parse newicks?
class NewickFormatter:
    def __init__(self, dist_matrix):
        self.dist_matrix = dist_matrix

    # in order for parse to work there needs to have been a tree object that is inserted into the class
    def parse_newick(self, newick):
        """takes a newick string and converts it into a simple binary tree with Biopythons phylo module"""
        tree = Phylo.read(StringIO(newick), "newick")
        return tree

        # place holder NJ algorithm formatting to Newick

    def generate_newick(self):
        # distance matrix and n x n dimensions, respectfully
        # dist_matrix = self.matrix
        tree = {}
        nj_nodes = [str(i) for i in range(len(self.dist_matrix))]
        while len(self.dist_matrix) != 2:
            self.dist_matrix, tree, nj_nodes = NeighborJoining._cluster_NJ(
                tree, nj_nodes
            )
            # merge remaining nodes in 2x2
        # perform merge on final 2 nodes
        node_a, node_b = nj_nodes[0], nj_nodes[1]
        dist = self.dist_matrix[0][1]
        limb_length = dist / 2

        # Make sure internal node is formatted correctly
        last_key = f"{node_a}<>{node_b}"
        tree[last_key] = {node_a: limb_length, node_b: limb_length}
        # clean up tree so nodes reflect <> notation
        final_tree = {}
        for key in tree:
            if key[0] == "[":
                final_tree[key[1 : len(key) - 1]] = tree[key]
            else:
                final_tree[key] = tree[key]
        return self._to_newick(final_tree)

    def _to_newick(self, tree):
        def recurse(node):
            if isinstance(tree[node], dict):
                children = []
                for child, dist in tree[node].items():
                    if child in tree:
                        children.append(f"{recurse(child)}:{dist}")
                    else:
                        children.append(f"{child}:{dist}")
                return f"({','.join(children)})"
            return node

        # Get the topmost node (root)
        root = list(tree.keys())[-1]
        return recurse(root) + ";"
