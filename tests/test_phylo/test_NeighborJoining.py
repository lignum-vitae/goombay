import unittest
import numpy as np
from goombay.phylo import NeighborJoining


class TestNeighborJoining(unittest.TestCase):
    def setUp(self):
        # Simple 3-taxa distance matrix
        self.dist_matrix_3 = np.array(
            [[0.0, 5.0, 9.0], [5.0, 0.0, 10.0], [9.0, 10.0, 0.0]], dtype=np.float32
        )

        # 4-taxa distance matrix
        self.dist_matrix_4 = np.array(
            [
                [0.0, 5.0, 9.0, 9.0],
                [5.0, 0.0, 10.0, 10.0],
                [9.0, 10.0, 0.0, 8.0],
                [9.0, 10.0, 8.0, 0.0],
            ],
            dtype=np.float32,
        )

    def test_total_row_distances(self):
        nj = NeighborJoining(self.dist_matrix_3)
        total = nj._total_row_distances()
        self.assertEqual(total, [14.0, 15.0, 19.0])

    def test_adjusted_distance(self):
        nj = NeighborJoining(self.dist_matrix_3)
        divergences = nj._total_row_distances()
        adj = nj._adjusted_distance(divergences)
        # Check symmetry and diagonal
        for i in range(3):
            self.assertEqual(adj[i][i], 0)
            for j in range(3):
                self.assertAlmostEqual(adj[i][j], adj[j][i])

    def test_pair_distance(self):
        nj = NeighborJoining(self.dist_matrix_3)
        pd = nj._pair_distance(0, 1)
        # Only one other node (2), so only one value
        self.assertEqual(len(pd), 1)
        self.assertAlmostEqual(pd[0], (9.0 + 10.0 - 5.0) / 2)

    def test_limb_length(self):
        nj = NeighborJoining(self.dist_matrix_3)
        divergences = nj._total_row_distances()
        dAZ, dBZ = nj._limb_length(0, 1, divergences)
        # Check that limb lengths sum to the distance between nodes
        self.assertAlmostEqual(dAZ + dBZ, nj.dist_matrix[0][1])

    def test_generate_newick_3taxa(self):
        nj = NeighborJoining(self.dist_matrix_3)
        newick = nj.generate_newick()
        # The output should be a valid Newick string ending with ';'
        self.assertTrue(newick.endswith(";"))
        self.assertIn("0", newick)
        self.assertIn("1", newick)
        self.assertIn("2", newick)
        self.assertEqual("(2:3.5,(0:2.0,1:3.0):3.5);", newick)

    def test_generate_newick_4taxa(self):
        nj = NeighborJoining(self.dist_matrix_4)
        newick = nj.generate_newick()
        self.assertTrue(newick.endswith(";"))
        for i in range(4):
            self.assertIn(str(i), newick)

    def test_generate_newick_2taxa(self):
        dist_matrix_2 = np.array([[0.0, 3.0], [3.0, 0.0]], dtype=np.float32)
        nj = NeighborJoining(dist_matrix_2)
        newick = nj.generate_newick()
        self.assertTrue(newick.endswith(";"))
        self.assertIn("0", newick)
        self.assertIn("1", newick)

    def test_to_newick_tree_format(self):
        # Test the _to_newick method with a simple tree dict
        nj = NeighborJoining(self.dist_matrix_3)
        tree = {
            "(1<>0)": {"0": 2.5, "1": 2.5},
            "((1<>0)<>2)": {"2": 4.5, "(1<>0)": 4.5},
        }
        newick = nj._to_newick(tree)
        self.assertTrue(newick.startswith("("))
        self.assertTrue(newick.endswith(";"))
        self.assertEqual(newick, "(2:4.5,(0:2.5,1:2.5):4.5);")


if __name__ == "__main__":
    unittest.main()
