########################################################################
##Aim of the script : use fastx to remove the adaptors and get QC plots##
##for 1 sample                                                         ##
#########################################################################

# module load fastx_toolkit

###########
##Library##
###########
import sys
import os, errno
import pandas as pd
import argparse
import subprocess
import gzip

#############
##Functions##
#############


def __create_necessary_output_directory(output_directory):
    out_path = os.path.join(output_directory, "ClippedFiles")
    try:
        os.mkdir(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    temp_path = os.path.join(output_directory, "TEMPFILES")
    try:
        os.mkdir(temp_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    QC_path = os.path.join(output_directory, "QC")
    try:
        os.mkdir(QC_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    stats_path = os.path.join(output_directory, "number_different_reads")
    try:
        os.mkdir(stats_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def __get_sample_information(sample_number):
    samples_informations = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/01.description_of_raw_data/locations_and_informations_of_fastq_files.tsv"
    )
    sample_information = samples_informations.iloc[sample_number]
    return sample_information


def __treat_input_variables():
    parser = argparse.ArgumentParser(
        description="Process the single integer expected in input"
    )
    parser.add_argument("row_number", nargs=1, type=int)
    try:
        row_number = parser.parse_args().row_number[0]
        return row_number
    except:
        print("Unrecognised input, this script is suppose to take a single integer")


def __create_study_files(df_row, data_output_directory):
    fastqgz_path = df_row.path_to_fastq

    fastq_path = os.path.join(
        data_output_directory,
        "TEMPFILES",
        df_row.batch_name + "_" + df_row.library_ID + ".fastq",
    )

    # First we unzip the corresponding files
    try:
        print("gunziping")
        print(fastqgz_path)
        with open(fastq_path, "w") as out_handle:
            p = subprocess.Popen(["gunzip", "-c", fastqgz_path], stdout=out_handle)
            p.communicate()
    except:
        raise

    # we output the QC files
    qc_path = os.path.join(
        data_output_directory,
        "QC",
        df_row.batch_name + "_" + df_row.library_ID + ".quality.txt",
    )
    try:
        print("computing QC stats")
        with open(qc_path, "w") as out_handle:
            p = subprocess.Popen(
                ["fastx_quality_stats", "-Q33", "-i", fastq_path], stdout=out_handle
            )
            p.communicate()
    except:
        raise

    # let's do the different clipping
    adaptor_sequence = "TGGAATTCTCGGGTGCCAAGGAACTCCAG"
    clipped_path = os.path.join(
        data_output_directory,
        "ClippedFiles",
        df_row.batch_name + "_" + df_row.library_ID + ".clipped.fastq",
    )
    try:
        print("clipping, keeping clipped")
        with open(clipped_path, "w") as out_handle:
            p = subprocess.Popen(
                [
                    "fastx_clipper",
                    "-a",
                    adaptor_sequence,
                    "-Q33",
                    "-l",
                    "0",
                    "-c",
                    "-n",
                    "-M",
                    "10",
                    "-i",
                    fastq_path,
                ],
                stdout=out_handle,
            )
            p.communicate()
            # -Q33 minimum quality asked
            # -l 0 we keep all sequence, no matter the length
            # -c Discard non-clipped sequences (i.e. - keep only sequences which contained the adapter)
            # -n keep sequences with unknown (N) nucleotides.
            # -M 10 require minimum adapter alignment length of 10
    except:
        raise

    clipped_adaptor_path = os.path.join(
        data_output_directory,
        "ClippedFiles",
        df_row.batch_name + "_" + df_row.library_ID + ".adaptoronly.fastq",
    )
    try:
        print("clipping, keeping adaptors")
        with open(clipped_adaptor_path, "w") as out_handle:
            p = subprocess.Popen(
                [
                    "fastx_clipper",
                    "-a",
                    adaptor_sequence,
                    "-Q33",
                    "-l",
                    "0",
                    "-k",
                    "-n",
                    "-M",
                    "10",
                    "-i",
                    fastq_path,
                ],
                stdout=out_handle,
            )
            p.communicate()
            # -k Report Adapter-Only sequences.
    except:
        raise

    unclipped_path = os.path.join(
        data_output_directory,
        "ClippedFiles",
        df_row.batch_name + "_" + df_row.library_ID + ".notclipped.fastq",
    )
    try:
        print("clipping, keeping unclipped")
        with open(unclipped_path, "w") as out_handle:
            p = subprocess.Popen(
                [
                    "fastx_clipper",
                    "-a",
                    adaptor_sequence,
                    "-Q33",
                    "-l",
                    "0",
                    "-C",
                    "-n",
                    "-M",
                    "10",
                    "-i",
                    fastq_path,
                ],
                stdout=out_handle,
            )
            p.communicate()
            # -C Discard clipped sequences (i.e. - keep only sequences which did not contained the adapter).
    except:
        raise


def __count_number_of_read_and_repartition_of_length_from_fastq(path_to_fastq):
    print("counting reads")
    results = {}
    total_number_of_read = 0
    with open(path_to_fastq, "r") as fastq_handle:
        count = 0
        for current_line in fastq_handle:
            if count == 1:
                RNA_read = current_line.strip()
                total_number_of_read += 1
                RNA_len = len(RNA_read)
                if RNA_len in results.keys():
                    results[RNA_len] += 1
                else:
                    results[RNA_len] = 1
            count = (count + 1) % 4

    length_of_read = []
    corresponding_number_of_reads = []
    for len_RNA, number_reads in results.items():
        length_of_read.append(len_RNA)
        corresponding_number_of_reads.append(number_reads)
    results = pd.DataFrame(
        {
            "length_of_read": length_of_read,
            "number_of_read": corresponding_number_of_reads,
        }
    )
    return results


def __count_number_of_read(df_row, data_output_directory):
    fastq_path = os.path.join(
        data_output_directory,
        "TEMPFILES",
        df_row.batch_name + "_" + df_row.library_ID + ".fastq",
    )
    clipped_path = os.path.join(
        data_output_directory,
        "ClippedFiles",
        df_row.batch_name + "_" + df_row.library_ID + ".clipped.fastq",
    )
    clipped_adaptor_path = os.path.join(
        data_output_directory,
        "ClippedFiles",
        df_row.batch_name + "_" + df_row.library_ID + ".adaptoronly.fastq",
    )
    unclipped_path = os.path.join(
        data_output_directory,
        "ClippedFiles",
        df_row.batch_name + "_" + df_row.library_ID + ".notclipped.fastq",
    )

    repartition_number_of_reads_full = __count_number_of_read_and_repartition_of_length_from_fastq(
        fastq_path
    )
    repartition_number_of_reads_clipped = __count_number_of_read_and_repartition_of_length_from_fastq(
        clipped_path
    )
    repartition_number_of_reads_clipped_adaptor = __count_number_of_read_and_repartition_of_length_from_fastq(
        clipped_adaptor_path
    )
    repartition_number_of_reads_unclipped = __count_number_of_read_and_repartition_of_length_from_fastq(
        unclipped_path
    )

    fastq_path = os.path.join(
        data_output_directory,
        "number_different_reads",
        df_row.batch_name + "_" + df_row.library_ID + ".all_reads.count.tsv",
    )
    clipped_path = os.path.join(
        data_output_directory,
        "number_different_reads",
        df_row.batch_name + "_" + df_row.library_ID + ".clipped.count.tsv",
    )
    clipped_adaptor_path = os.path.join(
        data_output_directory,
        "number_different_reads",
        df_row.batch_name + "_" + df_row.library_ID + ".adaptoronly.count.tsv",
    )
    unclipped_path = os.path.join(
        data_output_directory,
        "number_different_reads",
        df_row.batch_name + "_" + df_row.library_ID + ".notclipped.count.tsv",
    )

    repartition_number_of_reads_full.to_csv(fastq_path, sep="\t", index=False)
    repartition_number_of_reads_clipped.to_csv(clipped_path, sep="\t", index=False)
    repartition_number_of_reads_clipped_adaptor.to_csv(
        clipped_adaptor_path, sep="\t", index=False
    )
    repartition_number_of_reads_unclipped.to_csv(unclipped_path, sep="\t", index=False)


def __delete_unecessary_files(df_row, data_output_directory):
    fastq_path = os.path.join(
        data_output_directory,
        "TEMPFILES",
        df_row.batch_name + "_" + df_row.library_ID + ".fastq",
    )
    clipped_adaptor_path = os.path.join(
        data_output_directory,
        "ClippedFiles",
        df_row.batch_name + "_" + df_row.library_ID + ".adaptoronly.fastq",
    )
    unclipped_path = os.path.join(
        data_output_directory,
        "ClippedFiles",
        df_row.batch_name + "_" + df_row.library_ID + ".notclipped.fastq",
    )
    clipped_path = os.path.join(
        data_output_directory,
        "ClippedFiles",
        df_row.batch_name + "_" + df_row.library_ID + ".clipped.fastq",
    )

    os.remove(fastq_path)
    os.remove(clipped_adaptor_path)
    os.remove(unclipped_path)
    os.remove(clipped_path)


def __filter_clipped_on_size1826(df_row, data_output_directory):
    clipped_path = os.path.join(
        data_output_directory,
        "ClippedFiles",
        df_row.batch_name + "_" + df_row.library_ID + ".clipped.fastq",
    )
    clipped_filtered_path = os.path.join(
        data_output_directory,
        "ClippedFiles",
        df_row.batch_name + "_" + df_row.library_ID + ".clipped_filtered1826.fastq.gz",
    )
    group_of_four_lines = []
    with open(clipped_path, "r") as infile:
        with gzip.open(clipped_filtered_path, "w") as outfile:
            for line in infile:
                group_of_four_lines.append(line)
                if len(group_of_four_lines) == 4:
                    read_length = len(group_of_four_lines[1].strip())
                    if (read_length >= 18) and (read_length <= 26):
                        for l in group_of_four_lines:
                            outfile.write(l.encode("utf-8"))
                    group_of_four_lines = []


def main():
    # reading the inputs
    row_number = __treat_input_variables()

    # creating (if necessary) the output directory
    data_output_directory = (
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/"
    )
    __create_necessary_output_directory(data_output_directory)

    sample_informations = __get_sample_information(row_number)

    __create_study_files(sample_informations, data_output_directory)

    __count_number_of_read(sample_informations, data_output_directory)

    __count_number_of_read(sample_informations, data_output_directory)

    __filter_clipped_on_size1826(sample_informations, data_output_directory)

    __delete_unecessary_files(sample_informations, data_output_directory)


if __name__ == "__main__":
    main()
