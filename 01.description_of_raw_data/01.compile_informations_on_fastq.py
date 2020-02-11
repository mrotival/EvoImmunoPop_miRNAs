##############
##Librairies##
##############
import pandas as pd
from os import listdir
from os.path import isfile, join
import gzip

#############
##Functions##
#############


def get_folder_for_each_batch(dir):
    directories = [d for d in listdir(dir) if "ERC_Main_miRNAseq" in d]
    return directories


def treat_a_single_batch(dir_name, path_data, list_to_fill):
    # First we obtain the batch name and date
    batch_name = dir_name.split("_")[3]
    batch_date = dir_name.split("_")[4]

    # then we obtain the list of the fastq
    if batch_name == "3rdBatch":
        path_to_fastq_dirs = [
            join(path_data, dir_name, "05_FASTQ", "Sent_1"),
            join(path_data, dir_name, "05_FASTQ", "Sent_2"),
        ]
    elif batch_name == "1stBatch.pt2":
        path_to_fastq_dirs = [join(path_data, dir_name, "05_FASTQ_new48")]
    else:
        path_to_fastq_dirs = [join(path_data, dir_name, "05_FASTQ")]

    for path_to_fastq_dir in path_to_fastq_dirs:
        print(path_to_fastq_dir)
        fastq_files = [f for f in listdir(path_to_fastq_dir) if "fastq.gz" in f]

        for fastq_file in fastq_files:
            temp = fastq_file.split("-")
            individual = temp[0]
            condition = temp[1].split("_")[0]
            replicate = temp[1].split("_")[1].split(".")[0]
            library_ID = fastq_file.split(".")[0]
            # for each fastq we read the name of the machine

            with gzip.open(join(path_to_fastq_dir, fastq_file), "r") as fastq_handle:
                raw_fastq_info = next(fastq_handle)
                raw_fastq_info = raw_fastq_info.decode("utf-8").strip()
                sequencing_machine = ":".join(raw_fastq_info[1:].split(":")[:3])

            path_to_fastq = join(path_to_fastq_dir, fastq_file)
            list_to_fill.append(
                pd.DataFrame(
                    {
                        "individual": individual,
                        "condition": condition,
                        "library_number": replicate,
                        "library_ID": library_ID,
                        "batch_name": batch_name,
                        "batch_date": batch_date[4:],
                        "sequencing_machine": sequencing_machine,
                        "path_to_fastq": path_to_fastq,
                    },
                    index=[0],
                )
            )


def main():
    path_to_data_freeze = "/pasteur/projets/policy01/evo_immuno_pop/DATA_FREEZE/"
    directories = get_folder_for_each_batch(path_to_data_freeze)
    main_results = []
    for d in directories:
        treat_a_single_batch(d, path_to_data_freeze, main_results)
    main_results = pd.concat(main_results).reset_index()
    main_results.to_csv(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/01.description_of_raw_data/locations_and_informations_of_fastq_files.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
