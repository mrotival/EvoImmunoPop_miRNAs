##############
##Librairies##
##############
import pandas as pd

#############
##Functions##
#############


def get_informations_on_samples():
    samples_informations = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/01.description_of_raw_data/locations_and_informations_of_fastq_files.tsv"
    )
    return samples_informations


def main():
    samples_informations = get_informations_on_samples()
    grouped_samples_informations = (
        samples_informations.groupby(["batch_name", "sequencing_machine"])[
            "path_to_fastq"
        ]
        .nunique()
        .reset_index(name="count")
    )
    print(grouped_samples_informations)

    grouped_samples_informations.to_csv(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/01.description_of_raw_data/number_fastq_files_per_batch_and_sequencing_machine.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
