try:
    # external dependencies
    import numpy
    from numpy._typing import NDArray
except ImportError:
    raise ImportError("Numpy is not installed. Please pip install numpy to continue.")

try:
    from io import StringIO
    from Bio import Phylo
    from Bio.Phylo.BaseTree import Tree
except ImportError:
    raise ImportError(
        "Biopython is not installed. Please pip install biopython to continue"
    )


# takes a matrix from a clustering algorithm and outputs a newick tree, can also parse newicks?
class NewickFormatter:
    def __init__(self, dist_matrix: NDArray):
        self.dist_matrix = dist_matrix

    # in order for parse to work there needs to have been a tree object that is inserted into the class
    def parse_newick(self, newick: str) -> Tree:
        """takes a newick string and converts it into a simple binary tree with Biopythons phylo module"""
        tree = Phylo.read(StringIO(newick), "newick")
        return tree


class NeighborJoining:
    def __init__(self, dist_matrix: NDArray):
        self.dist_matrix = numpy.array(dist_matrix)

    # distance calculation for NJ
    # returns a list of total distances for forming Distance Matrix Prime
    def _total_row_distances(self):
        total_distances = []
        for i in range(len(self.dist_matrix)):
            row_sum = sum(self.dist_matrix[i])
            total_distances.append(row_sum)
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

    def _pair_distance(self, nodeI: int, nodeJ: int) -> list[numpy.float32]:
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

    # limb length is calculated slightly differently by taking the delta between
    # nodes A and B into consideration instead of divergences.
    # Calculate limb lengths for each leaf that is joined
    # return a tuple containing two values for each distance
    def _limb_length(
        self, nodeA: int, nodeB: int, divergences: list[numpy.float32]
    ) -> tuple[numpy.float32, numpy.float32]:
        n = len(self.dist_matrix)
        dAB = self.dist_matrix[nodeA][nodeB]
        divergenceA = divergences[nodeA]
        divergenceB = divergences[nodeB]
        deltaAB = (divergenceA - divergenceB) / (n - 2)
        # limb lengths
        dAZ = (dAB + deltaAB) / 2
        dBZ = (dAB - deltaAB) / 2

        return dAZ, dBZ

    def _to_newick(self, tree: dict[str, dict[str, numpy.float32]]) -> str:
        def recurse(node: str):
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

    def _cluster_NJ(self, tree: dict[str, dict[str, numpy.float32]], nodes: list[str]):
        mat_len = len(self.dist_matrix)

        if mat_len == 2:
            return self.dist_matrix, tree, nodes

        divergences = self._total_row_distances()
        adj_distance_matrix = self._adjusted_distance(divergences)
        min_val = float("inf")
        min_i, min_j = 0, 0

        # Find the pair with the minimum adjusted distance
        for i in range(mat_len):
            for j in range(mat_len):
                if i != j:
                    val = adj_distance_matrix[i][j]
                    if val < min_val:
                        min_val = val
                        min_i = i
                        min_j = j
        # Calculate limb lengths for the new node
        new_limbs = self._limb_length(min_i, min_j, divergences)
        new_limb_MI, new_limb_MJ = new_limbs[0], new_limbs[1]
        # Create new node label
        new_node = f"({nodes[min_j]}<>{nodes[min_i]})"
        # Add new node to tree
        tree[new_node] = {
            nodes[min_i]: new_limb_MI,
            nodes[min_j]: new_limb_MJ,
        }
        # Remove merged nodes and add new node
        nodes_to_remove = [nodes[min_i], nodes[min_j]]
        new_nodes = [n for n in nodes if n not in nodes_to_remove]
        new_nodes.append(new_node)
        # Calculate new distances for the new node to remaining nodes
        new_node_distances = self._pair_distance(min_i, min_j)
        # Build new distance matrix
        new_mat_len = mat_len - 1
        new_distance_matrix = numpy.zeros(
            (new_mat_len, new_mat_len), dtype=numpy.float32
        )
        # Fill in the new distance matrix
        idx = 0
        for i in range(mat_len):
            if i == min_i or i == min_j:
                continue
            jdx = 0
            for j in range(mat_len):
                if j == min_i or j == min_j:
                    continue
                new_distance_matrix[idx][jdx] = self.dist_matrix[i][j]
                jdx += 1
            idx += 1
        # Add distances for the new node
        for i in range(new_mat_len - 1):
            dist = new_node_distances[i]
            new_distance_matrix[i][new_mat_len - 1] = dist
            new_distance_matrix[new_mat_len - 1][i] = dist
        # Recursively cluster
        self.dist_matrix = new_distance_matrix
        return self._cluster_NJ(tree, new_nodes)

    def generate_newick(self) -> str:
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
