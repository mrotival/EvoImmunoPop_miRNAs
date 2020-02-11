import pandas as pd


def main():
    print("Hello World !")
    informations_on_fastq = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/filtered_and_combined_samples.tsv"
    )
    grouped_info = []
    for index, sample_to_treat in informations_on_fastq.iterrows():
        grouped_info.append(
            pd.read_table(
                "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_raw_counts/{}.miRNA_raw_counts.tsv".format(
                    sample_to_treat.library_ID
                )
            )
        )
    grouped_info = pd.concat(grouped_info)
    grouped_info.reset_index(drop=True, inplace=True)
    grouped_info.to_csv(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/miRNA_raw_counts.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
