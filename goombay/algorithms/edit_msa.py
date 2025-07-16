#built-in
from __future__ import annotations

#internal dependencies
from goombay.algorithms.editdistance import (
        Needleman_Wunsch, 
        Lowrance_Wagner, 
        Wagner_Fischer,
        Waterman_Smith_Beyer,
        Gotoh,
        Hirschberg,
        Jaro,
        Jaro_Winkler)

try:
    # external dependencies
    import numpy
    from numpy import float64
    from numpy._typing import NDArray
    #global packages serving as a placeholder for parsing newick strings - adahik
    from Bio import Phylo
    from io import StringIO
except ImportError:
    raise ImportError("Please pip install all dependencies from requirements.txt!")

def main():
    seq1 = "HOUSEOFCARDSFALLDOWN"
    seq2 = "HOUSECARDFALLDOWN"
    seq3 = "FALLDOWN"
    seq_list = [seq1, seq2, seq3]

    #print(feng_doolittle.align(seq1, seq2, seq3)) 
    #print(Feng_Doolittle.pairwise)
    #print(Feng_Doolittle.supported_pairwise_algs())
    output = Feng_Doolittle(seq_list)
    print(output())

class Feng_Doolittle():
    supported_pairwise = {
        "needleman_wunsch"    : Needleman_Wunsch,
        "jaro"                : Jaro,
        "jaro_winkler"        : Jaro_Winkler,
        "gotoh"               : Gotoh,
        "wagner_fischer"      : Wagner_Fischer,
        "waterman_smith_beyer": Waterman_Smith_Beyer,
        "hirschberg"          : Hirschberg,
        "lowrance_wagner"     : Lowrance_Wagner
    }

    abbreviations = {
        "nw" : "needleman_wunsch",
        "j"  : "jaro",
        "jw" : "jaro_winkler",
        "g"  : "gotoh",
        "wf" : "wagner_fischer",
        "wsb": "waterman_smith_beyer",
        "h"  : "hirschberg",
        "lw" : "lowrance_wagner"
    }
    def __init__(self, seqs, pairwise: str = "needleman_wunsch"):
        """Initialize Feng-Doolittle algorithm with chosen pairwise method"""
        # Get pairwise alignment algorithm
        if pairwise in self.supported_pairwise:
            self.pairwise = self.supported_pairwise[pairwise]()
        elif pairwise in self.abbreviations:
            self.pairwise = self.supported_pairwise[self.abbreviations[pairwise]]()
        else:
            raise ValueError(f"Unsupported pairwise alignment method: {pairwise}")
        self.seqs = seqs # list of sequences

    @classmethod
    def supported_pairwise_algs(cls):
        return list(cls.supported_pairwise)

    # helper functions
    def parse_newick(self, newick):
        """takes a newick string and converts it into a simple binary tree with Biopythons phylo module"""
        tree = Phylo.read(StringIO(newick), "newick")
        return tree
    
    def merge_profiles(self, profile1: list[str], profile2: list[str]) -> list[str]:
        # Pick first seq from each profile as representative (simplified for now)
        rep1 = profile1[0]
        rep2 = profile2[0]

        # Align the two representative sequences
        score_matrix, pointer = self.pairwise(rep1, rep2) # Andrew Dahik: Score matrix is here if one wanted to see it.
        aligned_rep1, aligned_rep2 = self.traceback(rep1, rep2, pointer)

        # Helper: apply gaps to all sequences in a profile
        def apply_gaps(profile, aligned_rep):
            gapped_profile = []
            for seq in profile:
                gapped_seq = []
                seq_i = 0
                for char in aligned_rep:
                    if char == "-":
                        gapped_seq.append('-')
                    else:
                        gapped_seq.append(seq[seq_i])
                        seq_i += 1
                gapped_profile.append(''.join(gapped_seq))
            return gapped_profile
        # Apply alignment gap pattern to all sequences
        aligned_profile1 = apply_gaps(profile1, aligned_rep1)
        aligned_profile2 = apply_gaps(profile2, aligned_rep2)

        return aligned_profile1 + aligned_profile2
    
    def traceback(self, query_seq: list[str], subject_seq: list[str], pointer) -> tuple[str, str]:
        aligned_qs = []
        aligned_ss = []
        i = len(query_seq) - 1
        j = len(subject_seq) - 1
        while i > 0 or j > 0:
            direction = pointer[i][j]
            if direction == 2: # diagnonal (match/mismatch)
                aligned_qs.append(query_seq[i])
                aligned_ss.append(subject_seq[j])
                i -= 1
                j -= 1
            elif direction == 3: # up (gap in subject)
                aligned_qs.append(query_seq[i])
                aligned_ss.append("-")
                i -= 1
            elif direction == 4: # left (gap in query)
                aligned_qs.append("-")
                aligned_ss.append(query_seq[j])
                j -= 1
            elif direction >= 5 and direction < 8:  # diagonal + up
                aligned_qs.append(query_seq[i])
                aligned_ss.append(subject_seq[j])
                i -= 1
                j -= 1
            else:
                break  # Should not happen

        # Append remaining character (the 0th one) if missed
        if i == 0 and j == 0:
            aligned_qs.append(query_seq[0])
            aligned_ss.append(subject_seq[0])
        elif i == 0:
            aligned_qs.append(query_seq[0])
            aligned_ss.append("-")
        elif j == 0:
            aligned_qs.append("-")
            aligned_ss.append(subject_seq[0])
        return ''.join(reversed(aligned_qs)), ''.join(reversed(aligned_ss))
    
    def __call__(self):
        """"""
        # This sets the unnormalized sequence distance
        seq_dist_matrix = numpy.zeros((len(self.seqs),len(self.seqs)), dtype=float64)
        profile_dict = {}
        aligned_map = {}
        #print(seq_dist_matrix)
        for i in range(len(self.seqs)):
            i_seq = self.seqs[i]
            profile_dict[str(i)] = [self.seqs[i]] #storing lists instead of strings
            for j in range(len(self.seqs)):
                if seq_dist_matrix[i][j] == 0:
                    if j != i:
                        j_seq = self.seqs[j]
                        alignment = self.pairwise
                        score_matrix, pointer = alignment(i_seq, j_seq)
                        aligned_qs, aligned_ss = self.traceback(i_seq, j_seq, pointer)
                        aligned_map[i_seq+"-"+j_seq] = (aligned_qs, aligned_ss)
                        print(aligned_map)
                        alignment_score = score_matrix[-1][-1]
                        seq_dist_matrix[i][j] = alignment_score
                        seq_dist_matrix[j][i] = alignment_score
                        #print(traceback)
        newick_tree = self.parse_newick(NeighborJoining(seq_dist_matrix).generateNJNewick())
        for clade in newick_tree.get_nonterminals(order='postorder'):
           
            left, right = clade.clades
            print(left, right)
            if left.name in profile_dict and right.name in profile_dict:
                print(f"Merging {left.name} and {right.name}")
                # Merge the aligned profiles
                left_profile = profile_dict[left.name]
                right_profile = profile_dict[right.name]
                merged_profile = self.merge_profiles(left_profile, right_profile) # these are consensus merges
                #print(merged_profile)
                # Assign unique name to internal node if needed
                if not clade.name:
                    clade.name = f"merged_{id(clade)}"
                #store the merged profile
                profile_dict[clade.name] = merged_profile

                # Optionally track who was merged
                #alignment_trace.append((left.name, right.name, clade.name))
        print(newick_tree)

        #print("Trace of merges:", alignment_trace)
        keys = list(profile_dict.keys())
        aligned_seqs = profile_dict[keys[-1]]
        print(aligned_seqs)
        rtn_str = ''
        for i in range(len(aligned_seqs)):
            if i != len(aligned_seqs) - 1:
                rtn_str = rtn_str + aligned_seqs[i] + "\n"
            else:
                rtn_str = rtn_str + aligned_seqs[i]
        return rtn_str
    
class NeighborJoining():
    def __init__(self, matrix):
        self.matrix = matrix
     # // distance calculation stuff for NJ
    #returns a list of total distances for forming Distance Matrix Prime    
    def _totalRowDistances(self):
        n = len(self.matrix)
        TotalDistances = []
        for i in range(n):
            TotalSum = 0
            for j in range(n):
                TotalSum += self.matrix[i][j]
            TotalDistances.append(TotalSum)
        return TotalDistances

    #adjustedDistanceFollows a different calculation instead
    def _adjustedDistance(self, divergences):
        adjMatrix = []
        L = len(self.matrix)
        for i in range(L):
            nodeRow = []
            adjMatrix.append(nodeRow)
            for j in range(L):
                if j == i:
                    nodeRow.append(0)
                else:
                    dIJ = self.matrix[i][j]
                    dIJ_prime = ((L-2)*dIJ) - (divergences[i] + divergences[j])
                    nodeRow.append(dIJ_prime)
        return adjMatrix
                    
    def _pairDistance(self, nodeI, nodeJ):
        #return new calculated distances 
        stored_values = []
        L = len(self.matrix)
        for k in range(L):
            if k != nodeI and k != nodeJ:
                #dMI/dMJ
                dM = (self.matrix[nodeI][k] + self.matrix[nodeJ][k] - self.matrix[nodeI][nodeJ])/2
                stored_values.append(dM)
                
        return stored_values

    #limb length is calculated slightly difference by taking the delta between
    #nodes A and B into consideration instead of divergences
    #calculate limb lengths for each leaf that is joined
    #return a tuple containing two values for each distance
    def _limbLength(self, nodeA, nodeB, divergences):

        n = len(self.matrix)
        dAB = self.matrix[nodeA][nodeB]
        divergenceA = divergences[nodeA]
        divergenceB = divergences[nodeB]
        deltaAB = (divergenceA - divergenceB) / (n-2)
        #limb lengths
        dAZ = (dAB + deltaAB)/2
        dBZ = (dAB - deltaAB)/2
        
        return dAZ, dBZ
    
        # // Cluster NJ Stuff
    def _toNewick(self, tree):
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

    def _clusterNJ(self, tree, nodes):
   
        L = len(self.matrix)
        if L == 2:
            return self.matrix
        divergences = self._totalRowDistances()
        ##print("total distances", divergences)
        adj_distance_matrix = self._adjustedDistance(divergences)
        minVal = float("inf")
        #looking for neighboring nodes,  based on minimum value, node i and node j
        #will be neighbors
        for i in range(L):
            for j in range(L):
                
                Val = adj_distance_matrix[i][j]
                if Val < minVal:
                    minVal = Val
                    #store node indices
                    I = i
                    J = j        
        #merge node I and node J together
        #grab distances of other nodes
        new_limbs = self._limbLength(I, J, divergences)
        #these limbs are for tree construction
        new_limb_MI, new_limb_MJ = new_limbs[0], new_limbs[1]
        #styling choice for tree
        tree[f"({nodes[J]}<>{nodes[I]})"] = {nodes[I]: new_limb_MI, nodes[J]: new_limb_MJ}
        nodes.insert(nodes.index(nodes[I]), f"({nodes[J]}<>{nodes[I]})")
        nodes.remove(str(nodes[I+1]))
        nodes.remove(str(nodes[J]))

        # calculate new distances for new node to remaining nodes
        # construct a new distance matrix
        # store distance matrice values for new node
        new_node_distances = self._pairDistance(I, J)

        # filler node that will be not used, only there to skip iteration
        new_distance_matrix = [[0 for _ in range(L-1)] for _ in range(L-1)]
        saved_matrices_values = []

        # k is the row_idx and m is the col_idx, change if desired
        # grabbing values that have not been merged
        for i in range(len(new_node_distances)):
            saved_matrices_values.append(new_node_distances[i])        
        for k in range(L):
            if k != I and k != J:
                for m in range(L-1):
                    if m != I and m != J:
                        if self.matrix[k][m] != 0:
                            saved_matrices_values.append(self.matrix[k][m])
        if L - 1 > 2:    
            # Fill in newDistanceMatrix based on conditions
            count = len(saved_matrices_values) - 1
            for i in range(L):
                
                for j in range(i+1, L):
                    if j - 1 == i:
                        new_distance_matrix[i][j-1] = 0
                        new_distance_matrix[j-1][i] = 0
                    else:
                        new_distance_matrix[i][j-1] = saved_matrices_values[j-count]
                        new_distance_matrix[j-1][i] = saved_matrices_values[j-count]
                count -= 1
                #now solve i for when j does not equal zero
        elif L - 1 == 2:
            new_distance_matrix[0][1] = saved_matrices_values[0]
            new_distance_matrix[1][0] = saved_matrices_values[0]
        #print(new_distance_matrix)
        return new_distance_matrix, tree, nodes
    # place holder NJ algorithm formatting to Newick
    def generateNJNewick(self):
        #distance matrix and n x n dimensions, respectfully
        dist_matrix = self.matrix
        tree = {}
        nj_nodes = [str(i) for i in range(len(dist_matrix))]
        while len(dist_matrix) != 2:
            dist_matrix, tree, nj_nodes = self._clusterNJ(tree, nj_nodes)
            #merge remaining nodes in 2x2
        # perform merge on final 2 nodes
        node_a, node_b = nj_nodes[0], nj_nodes[1]
        dist = dist_matrix[0][1]
        limb_length = dist / 2

        # Make sure internal node is formatted correctly
        last_key = f"{node_a}<>{node_b}"
        tree[last_key] = {node_a: limb_length, node_b: limb_length}
        #clean up tree so nodes reflect <> notation
        final_tree = {}
        for key in tree:
            if key[0] == "[":
                final_tree[key[1:len(key)-1]] = tree[key]
            else:
                final_tree[key] = tree[key]
        #print(final_tree)
        #newick_tree = toNewick(final_tree)
        #print(fitFM(final_tree, backup_nj_nodes))
        return self._toNewick(final_tree)

    
#feng_doolittle = Feng_Doolittle()

if __name__ == "__main__":
    main()
