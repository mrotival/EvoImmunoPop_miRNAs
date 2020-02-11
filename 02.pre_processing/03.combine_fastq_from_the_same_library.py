#############################################################################
##This script aims to go over the fastq files and combine those coming from##
##the same library                                                         ##
#############################################################################


############
##Packages##
############
import os
import subprocess
import pandas as pd

#############
##Functions##
#############


def main():
    # First we get the informations for each fastq
    samples_informations = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/01.description_of_raw_data/locations_and_informations_of_fastq_files.tsv"
    )

    # we check for each library ID wethr there is only 1 sample or several
    grouped = samples_informations.groupby("library_ID")
    print(samples_informations.head())
    combine_samples_information_list = []
    for name, group in grouped:
        temp = {}
        temp["individual"] = group.iloc[[0]].individual
        temp["condition"] = group.iloc[[0]].condition
        temp["library_number"] = group.iloc[[0]].library_number
        temp["library_ID"] = name
        combine_samples_information_list.append(pd.DataFrame(temp))

        out_file_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/filtered_and_combined_fastq/"
        out_file_path += "{}_filtered_combined.fastq".format(name)

        # print(name)
        # with open(out_file_path, "w") as out_file_handle :
        #     for index, row in group.iterrows() :
        #         infile_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/ClippedFiles/{}_{}.clipped_filtered1826.fastq.gz".format(row.batch_name, row.library_ID)
        #         p = subprocess.Popen(["gunzip", "-c", infile_path], stdout = out_file_handle)
        #         p.communicate()
        # p = subprocess.Popen(["gzip", out_file_path])
        # p.communicate()

    combine_samples_information = pd.concat(combine_samples_information_list)
    combine_samples_information.to_csv(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/filtered_and_combined_samples.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
#EUB136-4_1
