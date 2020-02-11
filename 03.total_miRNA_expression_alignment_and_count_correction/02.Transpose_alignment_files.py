######################################################################################################################################
##This script aims at getting the 5' 3' DNA sequence corresponding to the read, so to use the reverse complement in case of - strand##
##We also compute the depth after cmc and keep track of mismatches                                                                  ##
##All of this will be used for isomir computation                                                                                   ##
######################################################################################################################################

###########
##imports##
###########
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import pandas as pd
import argparse

#############
##functions##
#############
def __treat_input_variables(max_int):
    def check_within_bound(value):
        ivalue = int(value)
        if not ((ivalue >= 0) and (ivalue <= max_int)):
            raise argparse.ArgumentTypeError(
                "{} is not positive or above {}".format(value, max_int)
            )
        return ivalue

    print("treating input variables")
    parser = argparse.ArgumentParser(
        description="Process the single integer expected in input"
    )
    parser.add_argument("row_number", nargs=1, type=check_within_bound)
    row_number = parser.parse_args().row_number[0]
    return row_number


def transpose_file(infile_path, outfile_path):
    with open(infile_path, "r") as intersectFile, open(outfile_path, "w") as outFile:
        for line in intersectFile:
            line = line.split("\t")
            read_chromosome = line[0]
            read_start = int(line[1])
            read_end = int(line[2])
            read_sequence = Seq(line[3], generic_dna)
            read_strand = line[5]
            read_depth = int(line[6])
            alignment_weight = float(line[7])
            mismatches = line[8]
            miRNA_chromosome = line[9]
            miRNA_start = int(line[10])
            miRNA_end = int(line[11])
            miRNA_strand = line[14]
            miRNA_name = line[12]
            miRNA_name2 = line[15]
            miRNA_hairpin_name = line[16]

            weighted_depth = read_depth * alignment_weight
            final_seq = read_sequence
            mir_based_start = -300
            mir_based_end = -300

            try:
                assert read_strand == miRNA_strand
            except AssertionError:
                print("Bedtools error!")
                print(line)
                sys.exit(2)

            if read_strand == "-":
                mir_based_strat = miRNA_end - read_end
                mir_based_end = miRNA_start - read_start
                final_seq = str(read_sequence.reverse_complement())
            elif read_strand == "+":
                mir_based_strat = read_start - miRNA_start
                mir_based_end = read_end - miRNA_end
                final_seq = str(final_seq)
            else:
                print("WARNING NOT SUPPOSED TO HAPPEN")
                print(read_strand)

            new_line = "\t".join(
                [
                    read_chromosome,
                    str(mir_based_strat),
                    str(mir_based_end),
                    str(read_strand),
                    final_seq,
                    str(weighted_depth),
                    mismatches,
                    miRNA_name,
                    miRNA_name2,
                    miRNA_hairpin_name,
                ]
            )
            outFile.write(new_line)


def main():
    informations_on_fastq = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/filtered_and_combined_samples.tsv"
    )

    for index, sample_to_treat in informations_on_fastq.iterrows():
        ###################################################
        ##To be applied to intersection with maturemiRNA ##
        ###################################################
        print(index, sample_to_treat.library_ID)
        infile_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/alignment_on_mature_miRNA/bed_files/{}.maturemiRNA.intersect.bed".format(
            sample_to_treat.library_ID
        )
        outfile_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/alignment_on_mature_miRNA/after_transpose/{}.maturemiRNA.intersect.transposed.bed".format(
            sample_to_treat.library_ID
        )
        transpose_file(infile_path, outfile_path)


if __name__ == "__main__":
    main()
