##WARNING
##fastqc must be installed and in the path for this script to work
##If you are still on the TARS cluster, the command `module load fastqc/0.11.5` should do the trick

##############
##librairies##
##############
import pandas as pd
import subprocess
import os
import argparse

#############
##Functions##
#############


def get_informations_on_samples():
    samples_informations = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/01.description_of_raw_data/locations_and_informations_of_fastq_files.tsv"
    )
    return samples_informations


def compute_fastqc(row):
    fastq_gz_path = row.path_to_fastq
    out_file_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/scripts/01.description_of_raw_data/TEMP_FILES/{}_{}.fastq".format(
        row.batch_name, row.library_ID
    )
    print("gunziping")
    with open(out_file_path, "w") as out:
        subprocess.Popen(["gunzip", "-c", fastq_gz_path], stdout=out)
    print("computing QC")
    subprocess.call(
        [
            "fastqc",
            out_file_path,
            "-o",
            "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/01.description_of_raw_data/QC_results",
            "--extract",
            "-q",
        ]
    )
    os.remove(out_file_path)


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
    samples_informations = get_informations_on_samples()
    # row_number = __treat_input_variables(samples_informations.shape[0]-1)
    print(samples_informations.head())
    samples_informations = samples_informations[
        samples_informations["library_ID"].isin(
            ["EUB024-5_1", "EUB145-1_1", "EUB115-1_1", "EUB122-2_1"]
        )
    ]
    # row = samples_informations.iloc[row_number]

    for idx, row in samples_informations.iterrows():
        print(row.path_to_fastq)
        compute_fastqc(row)


if __name__ == "__main__":
    main()
