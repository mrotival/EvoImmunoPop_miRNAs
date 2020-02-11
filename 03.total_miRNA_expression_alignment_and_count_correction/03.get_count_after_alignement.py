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


def get_counts_for_one_sample(sample_infos):
    path_to_aligned_intersected_bed = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/alignment_on_mature_miRNA/bed_files/{}.maturemiRNA.intersect.bed".format(
        sample_infos.library_ID
    )
    alignments_info = pd.read_table(path_to_aligned_intersected_bed)
    alignments_info.columns = [
        "read_chromosome",
        "read_start",
        "read_end",
        "read_sequence",
        "V2",
        "read_strand",
        "read_depth",
        "alignement_weight",
        "misalignement",
        "miRNA_chromosome",
        "miRNA_start",
        "miRNA_end",
        "miRNA_name",
        "V3",
        "miRNA_strand",
        "miRNA_name_2",
        "miRNA_hairpin_name",
    ]
    to_iterate = alignments_info.groupby(
        ["miRNA_name", "miRNA_name_2", "miRNA_hairpin_name"]
    )
    individual_miRNA_count = []
    for miRNA_names, alignement_on_miRNA in to_iterate:
        # print(alignement_on_miRNA.head())
        number_counts = 0.0
        for idx, r in alignement_on_miRNA.iterrows():
            number_counts += r.read_depth + r.alignement_weight
        miRNA_count_one_sample = pd.DataFrame(
            {
                "individual": sample_infos.individual,
                "condition": sample_infos.condition,
                "library": sample_infos.library_number,
                "library_ID": sample_infos.library_ID,
                "miRNA_name": miRNA_names[0],
                "miRNA_name2": miRNA_names[1],
                "miRNA_hairpin": miRNA_names[2],
                "raw_counts": number_counts,
            },
            index=[0],
        )
        individual_miRNA_count.append(miRNA_count_one_sample)
    individual_miRNA_count = pd.concat(individual_miRNA_count)
    individual_miRNA_count.reset_index(drop=True, inplace=True)
    individual_miRNA_count.to_csv(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_raw_counts/{}.miRNA_raw_counts.tsv".format(
            sample_infos.library_ID
        ),
        sep="\t",
        index=False,
    )


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


def main():
    print("Hello World !")
    informations_on_fastq = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/filtered_and_combined_samples.tsv"
    )
    row_number = __treat_input_variables(informations_on_fastq.shape[0] - 1)
    sample_to_treat = informations_on_fastq.iloc[row_number]

    print(sample_to_treat.library_ID)
    get_counts_for_one_sample(sample_to_treat)


if __name__ == "__main__":
    main()
