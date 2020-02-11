#####################
##Aim of the script##
#####################
##We want to create the covariate files that will be use in
##several other scripts

###########
##Imports##
###########
import pandas as pd
import gzip

#############
##Functions##
#############


def treat_the_case_of_categorical_variables(samples_to_keep):
    print("treating categorical variables")
    categorical_variables_Katie_file = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Katie/ERC_Main_analyses/Covariates/AllCategoricalVariables.txt"
    )
    categorical_variables = categorical_variables_Katie_file.loc[samples_to_keep]
    categorical_variables.to_csv(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/covariates/all_categorical_variables.tsv",
        sep="\t",
    )


def get_gc(library_ID, all_fastq_informations):
    # Get the percent of GC from a library ID
    # We want to loop on all batch where this library ID appear
    temp = all_fastq_informations.copy()
    temp = temp[temp["library_ID"] == library_ID]
    GC_count = 0.0
    total_count = 0.0
    for idx, row in temp.iterrows():
        batch = row.batch_name
        gc_information = pd.read_table(
            "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/QC/{}_{}.quality.txt".format(
                batch, library_ID
            )
        )
        GC_count += gc_information["G_Count"].sum() + gc_information["C_Count"].sum()
        total_count += gc_information["Max_count"].sum()
    GC_percent = GC_count / total_count
    return GC_percent


def get_QC20(library_ID, all_fastq_informations):
    # Get the percent of GC from a library ID
    # We want to loop on all batch where this library ID appear
    temp = all_fastq_informations.copy()
    temp = temp[temp["library_ID"] == library_ID]
    QC20_count = 0.0
    qc_total_count = 0.0
    for idx, row in temp.iterrows():
        batch = row.batch_name
        QC_information_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/01.description_of_raw_data/QC_results/{}_{}_fastqc/fastqc_data.txt".format(
            batch, library_ID
        )
        with open(QC_information_path, "r") as QC_info:
            in_correct_bloc = False
            for l in QC_info:
                if ">>Per sequence quality scores" in l:
                    in_correct_bloc = True
                    continue
                if in_correct_bloc:
                    if "#Quality" in l:
                        continue
                    if ">>END_MODULE" in l:
                        in_correct_bloc = False
                        continue
                    else:
                        wline = l.strip().split("\t")
                        QC = int(wline[0])
                        count = float(wline[1])
                        qc_total_count += count
                        if QC > 20:
                            QC20_count += count
    return QC20_count / qc_total_count


def get_QC30(library_ID, all_fastq_informations):
    # Get the percent of GC from a library ID
    # We want to loop on all batch where this library ID appear
    temp = all_fastq_informations.copy()
    temp = temp[temp["library_ID"] == library_ID]
    QC30_count = 0.0
    qc_total_count = 0.0
    for idx, row in temp.iterrows():
        batch = row.batch_name
        QC_information_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/01.description_of_raw_data/QC_results/{}_{}_fastqc/fastqc_data.txt".format(
            batch, library_ID
        )
        with open(QC_information_path, "r") as QC_info:
            in_correct_bloc = False
            for l in QC_info:
                if ">>Per sequence quality scores" in l:
                    in_correct_bloc = True
                    continue
                if in_correct_bloc:
                    if "#Quality" in l:
                        continue
                    if ">>END_MODULE" in l:
                        in_correct_bloc = False
                        continue
                    else:
                        wline = l.strip().split("\t")
                        QC = int(wline[0])
                        count = float(wline[1])
                        qc_total_count += count
                        if QC > 30:
                            QC30_count += count
    return QC30_count / qc_total_count


def get_total_read(library_ID, all_fastq_informations):
    temp = all_fastq_informations.copy()
    temp = temp[temp["library_ID"] == library_ID]
    total_number_of_read = 0
    for idx, row in temp.iterrows():
        batch = row.batch_name
        total_number_read_file_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/number_different_reads/{}_{}.all_reads.count.tsv".format(
            batch, library_ID
        )
        with open(total_number_read_file_path, "r") as infile:
            next(infile)  # we skip the first line:
            for l in infile:
                total_number_of_read += int(l.strip().split()[1])
    return total_number_of_read


def get_clipped_read(library_ID, all_fastq_informations):
    temp = all_fastq_informations.copy()
    temp = temp[temp["library_ID"] == library_ID]
    clipped_number_of_read = 0
    for idx, row in temp.iterrows():
        batch = row.batch_name
        clipped_number_of_read_file_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/number_different_reads/{}_{}.clipped.count.tsv".format(
            batch, library_ID
        )
        with open(clipped_number_of_read_file_path, "r") as infile:
            next(infile)  # we skip the first line:
            for l in infile:
                clipped_number_of_read += int(l.strip().split()[1])
    return clipped_number_of_read


def get_adaptor_reads(library_ID, all_fastq_informations):
    temp = all_fastq_informations.copy()
    temp = temp[temp["library_ID"] == library_ID]
    clipped_number_of_adaptor = 0
    for idx, row in temp.iterrows():
        batch = row.batch_name
        clipped_number_of_adaptor_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/number_different_reads/{}_{}.adaptoronly.count.tsv".format(
            batch, library_ID
        )
        with open(clipped_number_of_adaptor_path, "r") as infile:
            next(infile)  # we skip the first line:
            for l in infile:
                clipped_number_of_adaptor += int(l.strip().split()[1])
    return clipped_number_of_adaptor


def get_notclipped_reads(library_ID, all_fastq_informations):
    temp = all_fastq_informations.copy()
    temp = temp[temp["library_ID"] == library_ID]
    number_of_not_clipped = 0
    for idx, row in temp.iterrows():
        batch = row.batch_name
        number_of_not_clipped_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/number_different_reads/{}_{}.notclipped.count.tsv".format(
            batch, library_ID
        )
        with open(number_of_not_clipped_path, "r") as infile:
            next(infile)  # we skip the first line:
            for l in infile:
                number_of_not_clipped += int(l.strip().split()[1])
    return number_of_not_clipped


def get_inserts1826(library_ID, all_fastq_informations):
    temp = all_fastq_informations.copy()
    temp = temp[temp["library_ID"] == library_ID]
    number_of_not_clipped = 0
    fastq_filtered_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/filtered_and_combined_fastq/{}_filtered_combined.fastq.gz".format(
        library_ID
    )
    number_line = 0.0
    with gzip.open(fastq_filtered_path, "r") as infile:
        for l in infile:
            number_line += 1
    return number_line / 4


def treat_the_case_of_continuous_variables(samples_to_keep, all_fastq_informations):
    print("treating continuous variables")
    print("getting possibles informations from Katie's files")
    continuous_variables_Katie_file = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Katie/ERC_Main_analyses/Covariates/AllContinuousVariables.txt"
    )
    column_to_keep = [
        "ng_ul",
        "ratio_260_280",
        "ratio_260_230",
        "rin",
        "quantity_sent",
        "cell_death",
        "total_axeq",
        "rin_axeq",
        "rrna_axeq",
    ]
    continuous_variable = continuous_variables_Katie_file.loc[
        samples_to_keep, column_to_keep
    ]

    print("iterating over samples to get other informations")
    continuous_variable_updated_list = []
    indexes = []
    test = 0
    for idx, row in continuous_variable.iterrows():
        test += 1
        print(idx, test)
        working_row = row.copy()
        working_row = working_row.append(
            pd.Series([get_gc(idx, all_fastq_informations)], index=["GC"])
        )
        working_row = working_row.append(
            pd.Series([get_QC20(idx, all_fastq_informations)], index=["QC20"])
        )
        working_row = working_row.append(
            pd.Series([get_QC30(idx, all_fastq_informations)], index=["QC30"])
        )
        working_row = working_row.append(
            pd.Series(
                [get_total_read(idx, all_fastq_informations)], index=["total_reads"]
            )
        )
        working_row = working_row.append(
            pd.Series(
                [get_clipped_read(idx, all_fastq_informations)], index=["clipped_reads"]
            )
        )
        working_row = working_row.append(
            pd.Series(
                [get_adaptor_reads(idx, all_fastq_informations)],
                index=["adaptor_reads"],
            )
        )
        working_row = working_row.append(
            pd.Series(
                [get_notclipped_reads(idx, all_fastq_informations)],
                index=["notclipped_reads"],
            )
        )
        # working_row = working_row.append(pd.Series([get_inserts1826(idx, all_fastq_informations)], index = ["inserts1826"]))

        indexes.append(idx)
        continuous_variable_updated_list.append(working_row)
    test = pd.DataFrame(continuous_variable_updated_list, index=indexes)
    test.to_csv(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/covariates/all_continuous_variables.tsv",
        sep="\t",
    )


def main():
    print("Hello, World !")
    print("getting informations on original fastq and fused fastq")
    # We need all fastq because it contains informations on batch
    all_fastq_informations = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/01.description_of_raw_data/locations_and_informations_of_fastq_files.tsv"
    )
    filtered_and_combined_fastq = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/filtered_and_combined_samples.tsv"
    )
    samples_to_keep = list(filtered_and_combined_fastq.library_ID.unique())
    samples_to_keep.sort()

    treat_the_case_of_categorical_variables(samples_to_keep)
    treat_the_case_of_continuous_variables(samples_to_keep, all_fastq_informations)


if __name__ == "__main__":
    main()
