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
    print(samples_informations.head())
    samples_informations["population"] = samples_informations.individual.apply(
        lambda x: x[:3]
    )
    samples_informations["individual_number"] = samples_informations.individual.apply(
        lambda x: x[3:]
    )
    samples_informations = samples_informations[
        ["population", "individual_number", "condition"]
    ]
    samples_informations = samples_informations.drop_duplicates()
    number_individuals_with_information = (
        samples_informations.groupby(["population", "condition"])["individual_number"]
        .nunique()
        .reset_index(name="count")
    )
    number_individuals_with_information.to_csv(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/01.description_of_raw_data/number_individuals_with_at_least_one_library.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
