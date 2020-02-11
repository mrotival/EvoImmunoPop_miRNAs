#####################################
##To read before runing this script##
#####################################
# bowtie must be in your path
# module load samtools/1.9
# module load bowtie/1.1.1
# module load Python/2.7.8
# module load bedtools/2.25.0
###########
##Imports##
###########
import pandas as pd
import argparse
import os
import gzip
import subprocess

#############
##Functions##
#############
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


def __create_necessary_folders(TEMP_folder):
    try:
        os.makedirs(TEMP_folder)
    except OSError:
        pass


def create_simplified_fasta(fastqgz_path, TEMP_folder, library_ID):
    """this function transform the fastq into a simplified fasta containing each version of the
    read exactly once"""

    print("writing simplified fasta")
    # First we get the counts
    count_of_each_read = {}
    with gzip.open(fastqgz_path, "r") as infile:
        group_of_four_lines = []
        for line in infile:
            group_of_four_lines.append(line)
            if len(group_of_four_lines) == 4:
                read_sequence = group_of_four_lines[1].decode("utf-8").strip()
                if read_sequence in count_of_each_read.keys():
                    count_of_each_read[read_sequence] += 1
                else:
                    count_of_each_read[read_sequence] = 1
                group_of_four_lines = []

    # Then we write the pseudo fasta
    index = 0
    simplified_fasta_path = os.path.join(
        TEMP_folder, "{}.trimmed_filtered.simplified.fasta".format(library_ID)
    )
    with open(simplified_fasta_path, "w") as simplified_fasta:
        for read_sequence, count in count_of_each_read.items():
            index += 1
            simplified_fasta.write(">Read{}:{}\n".format(index, count))
            simplified_fasta.write("{}\n".format(read_sequence))

    return simplified_fasta_path


def alignement_summary(path_to_sam_alignement, path_to_summary_file, library_ID):
    aligned_reads = {}
    unaligned_reads = {}

    with open(path_to_summary_file, "w") as summary_file:
        with open(path_to_sam_alignement, "r") as sam_file:
            for line in sam_file:
                if not (line[0] == "@"):  # we skip the header
                    wline = line.strip().split("\t")
                    read_name, read_count = wline[0].split(":")
                    read_count = int(read_count)
                    if wline[1] == "4":  # The read is not aligned
                        unaligned_reads[read_name] = read_count
                    else:  # The read is aligned at least once
                        aligned_reads[read_name] = read_count

            to_output = {
                "aligned_reads": 0,
                "aligned_sequences": 0,
                "unaligned_reads": 0,
                "unaligned_sequences": 0,
                "library_ID": library_ID,
            }

            for key, value in aligned_reads.items():
                to_output["aligned_reads"] += value
                to_output["aligned_sequences"] += 1

            for key, value in unaligned_reads.items():
                to_output["unaligned_reads"] += value
                to_output["unaligned_sequences"] += 1

            to_output = pd.DataFrame(to_output, index=[0])
            to_output.to_csv(path_to_summary_file, sep="\t", index=False)


def align_on_hg19(TEMP_folder, library_ID, bowtie_index, fasta_path):
    print("aligning on HG19 with bowtie")
    alignement_file_path = os.path.join(
        TEMP_folder, "{}.aligned_on_hg19.sam".format(library_ID)
    )
    alignement_summary_file_path = os.path.join(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/alignment_on_hg19_summary",
        "{}.alignement_on_hg_19_summary.tsv".format(library_ID),
    )
    print(
        [
            "bowtie",
            bowtie_index,
            "-f",
            "-a",
            "--best",
            "--strata",
            "-v",
            "2",
            "-m",
            "50",
            "--sam",
            fasta_path,
            alignement_file_path,
        ]
    )
    p = subprocess.Popen(
        [
            "bowtie",
            bowtie_index,
            "-f",
            "-a",
            "--best",
            "--strata",
            "-v",
            "2",
            "-m",
            "50",
            "--sam",
            fasta_path,
            alignement_file_path,
        ]
    )
    p.communicate()
    ################################################################################################
    ##	-v 2  No more than 2 mismatches                                                           ##
    ## 	-a we output all valid Alignments that valids the following options                       ##
    ##	-a --best same output but in the order best to worst                                      ##
    ##	-a --best --strata report only alignement in the best stratum (least number of mismatch)  ##
    ##	-m 50 do not report the reads that have been aligned to more than 50 places               ##
    ##  -f the input is fasta files                                                               ##
    ##	--sam prin alignment in SAM format                                                        ##
    ################################################################################################
    alignement_summary(alignement_file_path, alignement_summary_file_path, library_ID)
    return alignement_file_path


def format_alignement_for_cmc(TEMP_folder, library_ID, sam_file):
    align_almost_ready_for_cmc_path = os.path.join(
        TEMP_folder, "{}.cmc_almost_ready_alignement".format(library_ID)
    )
    with open(sam_file, "r") as infile:
        with open(align_almost_ready_for_cmc_path, "w") as outfile:
            for line in infile:
                if line[0] != "@":
                    wline = line.split("\t")
                    read_id, number_read = wline[0].split(":")
                    wline.insert(11, "XC:i:" + str(number_read))
                    if wline[5] != "*":
                        new_line = "\t".join(wline)
                        outfile.write(new_line)
                else:
                    outfile.write(line)

    return align_almost_ready_for_cmc_path


def pretreatment_suggested_by_cmc_manual(
    TEMP_folder, library_ID, align_almost_ready_for_cmc_path, fasta_human_genome
):
    print("pretreatment cmc")
    align_ready_for_cmc_path = os.path.join(
        TEMP_folder, "{}.cmc_ready_alignement.sam".format(library_ID)
    )
    align_sorted = os.path.join(
        TEMP_folder, "{}.cmc_ready_alignement.sorted.bam".format(library_ID)
    )

    p = subprocess.Popen(
        [
            "samtools",
            "import",
            fasta_human_genome,
            align_almost_ready_for_cmc_path,
            align_ready_for_cmc_path,
        ]
    )
    p.communicate()

    with open(align_sorted, "w") as outfile:
        p = subprocess.Popen(
            ["samtools", "sort", align_ready_for_cmc_path], stdout=outfile
        )
        p.communicate()

    with open(align_ready_for_cmc_path, "w") as outfile:
        p = subprocess.Popen(
            ["samtools", "fillmd", align_sorted, fasta_human_genome], stdout=outfile
        )
        p.communicate()

    return (align_ready_for_cmc_path, align_sorted)


def run_cmc(TEMP_folder, library_ID, aligned_sam_path):
    print("running cmc")
    cmc_path = "/pasteur/projets/policy01/evo_immuno_pop/Martin/Ressources/Programs/cmc-2010.03.15/cmc.py"
    weighted_alignement_path = os.path.join(
        TEMP_folder, "{}.aligned_on_hg_19.weighted.sam".format(library_ID)
    )
    p = subprocess.Popen(
        [
            "python2.7",
            cmc_path,
            "-i",
            aligned_sam_path,
            "-o",
            weighted_alignement_path,
            "-n",
            "10",
        ]
    )
    p.communicate()
    return weighted_alignement_path


def sam_to_bed(file_path):
    print("sam_to_bed")
    outfile_path = file_path[:-3] + "bed"
    with open(file_path, "r") as sam_file:
        with open(outfile_path, "w") as bed_file:
            for line in sam_file:
                if line[0] == "@":
                    continue
                wline = line.strip().split("\t")
                chromosome = "chr" + wline[2]
                if wline[1] == "16":
                    strand = "-"
                elif wline[1] == "0":
                    strand = "+"
                else:
                    print("not supposed to happen {}".format(line[1]))
                start = int(wline[3])
                end = start + len(wline[9]) - 1
                name = wline[9]  # it's the sequence
                ##We also need to keep track of the read depth and the weight of this alignment
                readDepth = wline[0].split(":")[1]
                weight = wline[15].split(":")[-1]
                ##And the nucleotide changes compared to reference
                mismatchingPositions = wline[13]
                new_line = (
                    "\t".join(
                        [
                            chromosome,
                            str(start),
                            str(end),
                            name,
                            ".",
                            strand,
                            readDepth,
                            weight,
                            mismatchingPositions,
                        ]
                    )
                    + "\n"
                )
                bed_file.write(new_line)
    return outfile_path


def intersection_with_mature_read(input_bed_file, output_bed_file):
    print("intersection with mature read")
    mirbase20_mature_path = "/pasteur/projets/policy01/evo_immuno_pop/ERCPilot_SharedDBs/mirbase20/miRNA_mature_coordinates_strandinfo.bed"
    with open(output_bed_file, "w") as output_bed:
        p = subprocess.Popen(
            [
                "bedtools",
                "intersect",
                "-a",
                input_bed_file,
                "-b",
                mirbase20_mature_path,
                "-s",
                "-f",
                "0.75",
                "-wb",
                "-wa",
            ],
            stdout=output_bed,
        )
        p.communicate()


def intersection_with_hairpin(input_bed_file, output_bed_file):
    print("intersection with hairpin")
    mirbase20_hairpin_path = "/pasteur/projets/policy01/evo_immuno_pop/ERCPilot_SharedDBs/mirbase20/miRNA_hairpin_coordinates_strandinfo.bed"
    with open(output_bed_file, "w") as output_bed:
        p = subprocess.Popen(
            [
                "bedtools",
                "intersect",
                "-a",
                input_bed_file,
                "-b",
                mirbase20_hairpin_path,
                "-s",
                "-f",
                "0.75",
                "-wb",
                "-wa",
            ],
            stdout=output_bed,
        )
        p.communicate()


def main():
    print("Starting alignment script")

    ##Setup before starting (variables definitions, creation of necessary folder, check which individual to use ...)
    fasta_human_genome = "/pasteur/projets/policy01/evo_immuno_pop/ERCPilot_SharedDBs/HumanGenomeGRCh37/GRCh37_dna_sm_ref.fa"
    bowtie_index_human_genome = "/pasteur/projets/policy01/evo_immuno_pop/ERCPilot_SharedDBs/HumanGenomeGRCh37/BowtieIndex/GRCh37_dna_sm_ref"
    informations_on_fastq = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/filtered_and_combined_samples.tsv"
    )
    row_number = __treat_input_variables(informations_on_fastq.shape[0] - 1)
    sample_to_treat = informations_on_fastq.iloc[row_number]
    fastqgz_location = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/02.pre_processing/filtered_and_combined_fastq/{}_filtered_combined.fastq.gz".format(
        sample_to_treat.library_ID
    )

    TEMP_folder = "/pasteur/scratch/mrotival/miRNA/temp_alignement_files"
    __create_necessary_folders(TEMP_folder)

    ########################
    ##Alignement on crgh37##
    ########################
    simplified_fasta_path = create_simplified_fasta(
        fastqgz_location, TEMP_folder, sample_to_treat.library_ID
    )
    align_hg19_sam_path = align_on_hg19(
        TEMP_folder,
        sample_to_treat.library_ID,
        bowtie_index_human_genome,
        simplified_fasta_path,
    )
    ###############################
    ##Use Cross Mapper Correction##
    ###############################

    align_almost_ready_for_cmc_path = format_alignement_for_cmc(
        TEMP_folder, sample_to_treat.library_ID, align_hg19_sam_path
    )
    align_ready_for_cmc_path, temp_must_be_deleted = pretreatment_suggested_by_cmc_manual(
        TEMP_folder,
        sample_to_treat.library_ID,
        align_almost_ready_for_cmc_path,
        fasta_human_genome,
    )
    weighted_alignement_sam = run_cmc(
        TEMP_folder, sample_to_treat.library_ID, align_ready_for_cmc_path
    )

    ############################################
    ##We check the intersection with mirbase20##
    ############################################
    alignement_weighted_bed = sam_to_bed(weighted_alignement_sam)

    bed_mature_miRNA_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/alignment_on_mature_miRNA/bed_files/{}.maturemiRNA.intersect.bed".format(
        sample_to_treat.library_ID
    )
    intersection_with_mature_read(alignement_weighted_bed, bed_mature_miRNA_path)

    bed_hairpin_miRNA_path = "/pasteur/projets/policy01/evo_immuno_pop/Maxime/miRNA_V2/data/03.total_miRNA_expression_alignment_and_count_correction/alignment_on_hairpin/bed_files/{}.hairpin.intersect.bed".format(
        sample_to_treat.library_ID
    )
    intersection_with_hairpin(alignement_weighted_bed, bed_hairpin_miRNA_path)

    # we delete all files that must be deleted
    files_to_delete = [
        simplified_fasta_path,
        align_hg19_sam_path,
        align_almost_ready_for_cmc_path,
        align_ready_for_cmc_path,
        weighted_alignement_sam,
        # alignement_weighted_bed, #We keep it in case of wanting to check for cross mapping
        temp_must_be_deleted,
    ]
    for f in files_to_delete:
        try:
            os.remove(f)
        except:
            pass

    print("alignment script ended succesfully")


if __name__ == "__main__":
    main()
