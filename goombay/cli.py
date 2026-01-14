# this will be the package call
from goombay.align.edit_msa import FengDoolittle, TCoffee


def main():
    # FengDoolittle
    # seq1 = "HOUSEOFCARDSFALLDOWN"
    # seq2 = "HOUSECARDFALLDOWN"
    # seq3 = "FALLDOWN"
    # seq_list = [seq1, seq2, seq3]

    # print(FengDoolittle.align(seq_list))
    # print(FengDoolittle.supported_pairwise_algs())
    # print(FengDoolittle.supported_clustering_algs())
    # fd_gotoh = FengDoolittle(pairwise="gotoh")

    # seqs from rna.informatik https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Notredame-Higgins-Heringa
    # seq_list = ["GARFIELD_THE_LAST_FAT_CAT", "GARFIELD_THE_FAST_CAT", "GARFIELD_THE_VERY_FAST_CAT", "THE_FAT_CAT"]

    # seqs from https://backofenlab.github.io/BioinformaticsII-pages/exercise-sheet-3.html
    seq_list = ["CACCGG", "ACCAAG", "AACACC"]
    # print(TCoffee.supported_pairwise_algs())
    # print(TCoffee.supported_clustering_algs())
    # print(t_coffee(seq_list))
    # fd_gotoh = TCoffee(pairwise="gotoh")
    print("=========FengDoolittle=========\n")
    feng_doo = FengDoolittle()
    feng_doo_list = feng_doo.get_matrix(seq_list)
    print(feng_doo_list)
    feng_doo_align = feng_doo.align(seq_list)
    print(feng_doo_align)
    print("\n=========TCoffee=========\n")
    t_coffee = TCoffee()
    t_coffee_list = t_coffee.get_matrix(seq_list)
    print(t_coffee_list)
    t_coffe_align = t_coffee.align(seq_list)
    print(t_coffe_align)


if __name__ == "__main__":
    main()
