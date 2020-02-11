######################################################################################
##Aim of the script : the same mutation tag (MD:Z:xxx) only refers to the REF genome ##
##However, if we want not to keep track of the entire transcript, we will also need ##
##to keep track of the transcript state                                             ##
######################################################################################
###########
##Imports##
###########
import pandas as pd
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein



#############
##Functions##
#############
def __treat_input_variables(max_int):
    """treat the variable given to the script
    This script should take only one int as a variable
    and this int shouldn't exceed the "max_int" value."""

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


def __create_new_tag(mir_seq, mir_strand, SAM_tag):
    """Does the shift in tag.
    input tags are in the format MD:Z:18G1T0 (for + strand)
    Wanted_tag will be in this format : 2;18A->G;20C->T;
    meaning there are 2 substitutions, one at the 18th nucleotide of the miRNA (not the REF)
    were the REF was *on the corresponding strand* A but became G """

    SAM_tag_interesting_part = SAM_tag.split(":")[2]
    if SAM_tag_interesting_part.isnumeric():
        # Mean there is no substitution
        return "0;"
    else:
        # First we create substrings in the type "20A"
        substring_list = []
        last_character_was_alpha = False
        current_substring = ""
        for chara in SAM_tag_interesting_part:
            if not last_character_was_alpha:
                current_substring += chara
                last_character_was_alpha = chara.isalpha()
            else:
                if chara.isalpha():
                    current_substring += chara
                else:
                    substring_list.append(current_substring)
                    current_substring = chara
                    last_character_was_alpha = False
        # Now that we have the substring list, we can more easily create the corresponding tag
        number_substitutions = len(substring_list)
        if mir_strand == "+":
            seq_positive_strand = mir_seq
        else:
            seq_positive_strand = mir_seq.reverse_complement()

        # Now we create the real thing
        new_tag = str(number_substitutions) + ";"
        if mir_strand == "+":
            seq_position = 0
            for substitution in substring_list:
                # hypothesis, only SNPs
                seq_position += int(substitution[:-1]) # position being changed on the + strand, 0-based
                ref_nucleotide = substitution[-1]
                new_nucleotide = seq_positive_strand.tostring()[seq_position]
                new_tag += (
                    str(seq_position + 1) + ref_nucleotide + "->" + new_nucleotide + ";"
                )
                seq_position += 1 # position being changed on the + strand, 1-based

        if mir_strand == "-":
            seq_position_positive_strand = 0
            tag_parts = []
            for substitution in substring_list:
                # hypothesis, only SNPs
                seq_position_positive_strand += int(substitution[:-1]) # position being changed on the + strand, 0-based
                ref_nucleotide_positive_strand = substitution[-1]
                new_nucleotide_positive_strand = str(seq_positive_strand)[
                    seq_position_positive_strand
                ]
                ref_nucleotide = str(
                    Seq(ref_nucleotide_positive_strand).reverse_complement()
                )
                new_nucleotide = str(
                    Seq(new_nucleotide_positive_strand).reverse_complement()
                )
                seq_position = ( 
                    len(str(seq_positive_strand)) - seq_position_positive_strand # position being changed on the - strand, 1-based
                )
                tag_parts.append(
                    str(seq_position) + ref_nucleotide + "->" + new_nucleotide + ";"
                )
                seq_position_positive_strand += 1 # position being changed on the + strand, 1-based
            tag_parts = tag_parts[::-1]
            for t in tag_parts:
                new_tag += t
        return new_tag


def treat_one_sample(sample_infos):
    """This function take one sample information on alignement and create a
    new column based on the strand, the subsitution and the sequence column.
    This new column will contain all the informations on the substitutions so that
    there is no need to know th sequence to classify isomirs"""

    input_file = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/alignment_on_mature_miRNA/after_transpose/{}.maturemiRNA.intersect.transposed.bed".format(
        sample_infos.library_ID
    )
    output_file = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03b.isomirs_alignment/alignement/after_transpose_and_substitutions/{}.maturemiRNA.intersect.transposed.substituted.bed".format(
        sample_infos.library_ID
    )

    with open(input_file, "r") as original_file, open(output_file, "w") as out_file:
        for line in original_file:
            wLine = line.strip().split()
            # print(wLine)
            strand = wLine[3]
            susbtitution_tag = wLine[6]
            sequence = Seq(wLine[4], generic_dna)
            new_tag = __create_new_tag(sequence, strand, susbtitution_tag)
            wLine[6] = new_tag
            out_file.write("\t".join(wLine) + "\n")


def main():
    informations_on_fastq = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/filtered_and_combined_samples.tsv"
    )

    for row_number in range(informations_on_fastq.shape[0]) :
        sample_to_treat = informations_on_fastq.iloc[row_number]
        print(sample_to_treat.library_ID)
        treat_one_sample(sample_to_treat)


if __name__ == "__main__":
    print("Hello World !")
    main()
    print("Script succesfully ended")
