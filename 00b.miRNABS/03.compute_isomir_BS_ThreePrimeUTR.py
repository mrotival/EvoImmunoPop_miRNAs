#############################################################################################################
##Aim of this script : compute the binding sites of the main 653 expressed miRNA (including those on X chr)##
##Onto the protein coding genes (not necessrily in 5'UTR)                                                  ##
#############################################################################################################

###########
##Imports##
###########

import pandas as pd
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
from subprocess import run
from tqdm import tqdm

#############
##Functions##
#############
def get_list_of_expressed_miRNAs():
    expression_matrix = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Martin/miRNA/03.alignment_and_count_correction/data/miRNA_counts.RPM.GC_Batch_corrected.tsv"
    )
    expressed_miRNAs = expression_matrix.index.values
    return expressed_miRNAs


def create_fasta_file_for_isomir(isomir_fasta_outfile):
    information_on_isomirs = pd.read_csv(
        "/pasteur/projets/policy01/evo_immuno_pop/Martin/miRNA/09.isomirs/data/isomiR_counts.RPM.GC_Batch_corrected.aggregated.tsv",
        sep="\t",
    )
    information_on_isomirs = information_on_isomirs[
        ["chromosome", "start", "end", "strand", "mir1", "mir2", "hairpin", "isomir_ID"]
    ]
    miRNA_coordinates = pd.read_csv(
        "/pasteur/projets/policy01/evo_immuno_pop/ERCPilot_SharedDBs/mirbase20/miRNA_mature_coordinates_strandinfo.bed",
        sep="\t",
    )
    miRNA_coordinates.columns = [
        "chromosome",
        "miRNA_start",
        "miRNA_end",
        "mir1",
        "v5",
        "strand",
        "mir2",
        "hairpin",
    ]

    information_on_isomir_with_known_coordinates = []
    for miRNA, isomir_of_miRNA in tqdm(information_on_isomirs.groupby(["mir1"])):
        isomir_info_expanded = isomir_of_miRNA.copy()
        temp = miRNA_coordinates[miRNA_coordinates.mir1 == miRNA]
        if temp.shape[0] == 1:
            isomir_info_expanded["miRNA_start"] = temp["miRNA_start"].values[0]
            isomir_info_expanded["miRNA_end"] = temp["miRNA_end"].values[0]
            isomir_info_expanded["chromosome"] = temp["chromosome"].values[0]
            information_on_isomir_with_known_coordinates.append(isomir_info_expanded)

    information_on_isomir_with_known_coordinates = pd.concat(
        information_on_isomir_with_known_coordinates
    )
    information_on_isomir_with_known_coordinates.reset_index(inplace=True, drop=True)

    human_genome = get_human_reference_genome()
    all_isomirs_record = []
    for row in tqdm(information_on_isomir_with_known_coordinates.itertuples()):
        if row.strand == "+":
            isomir_seq = (
                human_genome[row.chromosome[3:]]
                .seq[((row.miRNA_start - 1) + row.start) : (row.miRNA_end + row.end)]
                .transcribe()
            )
        else:
            isomir_seq = (
                human_genome[row.chromosome[3:]]
                .seq[((row.miRNA_start - 1) - row.end) : (row.miRNA_end - row.start)]
                .reverse_complement()
                .transcribe()
            )

        isomir_record = SeqRecord(
            isomir_seq,
            id="_".join([row.mir1, str(row.start), str(row.end), row.strand]),
            description="isomir_{}_{}_{}".format(
                row.mir1, row.start, row.end, row.strand
            ),
        )
        all_isomirs_record.append(isomir_record)
    with open(isomir_fasta_outfile, "w") as output_handle:
        SeqIO.write(all_isomirs_record, output_handle, "fasta")


def get_human_reference_genome():
    fasta_human_genome = "/pasteur/projets/policy01/evo_immuno_pop/ERCPilot_SharedDBs/HumanGenomeGRCh37/GRCh37_dna_sm_ref.fa"
    record_human = SeqIO.to_dict(SeqIO.parse(fasta_human_genome, "fasta"))
    return record_human


def get_RNA_of_each_transcript_in_fasta(output_path):
    ##First we read the informations we have about transcripts
    transcripts_info = pd.read_table(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/Evo_Immuno_pop_data/ExonCoordinates_hg37_ens70.txt",
        dtype="str",
    )
    transcripts_info["Chromosome Name"]
    transcripts_info.sort_values(by=["Chromosome Name"], inplace=True)
    # We only want to keep lines where chromosome is in the "chromosome_to_keep" list
    chromosome_to_keep = [str(k) for k in range(1, 23)]
    chromosome_to_keep.append("X")
    transcripts_info = transcripts_info[
        transcripts_info["Chromosome Name"].isin(chromosome_to_keep)
    ]
    transcripts_info = transcripts_info[~transcripts_info["3' UTR Start"].isnull()]
    transcripts_info = transcripts_info[~transcripts_info["3' UTR End"].isnull()]
    transcripts_info = transcripts_info[
        transcripts_info["Gene Biotype"] == "protein_coding"
    ]

    splitted_transcripts_info = transcripts_info.groupby(
        ["Ensembl Gene ID", "Ensembl Transcript ID", "Chromosome Name", "Strand"]
    )
    human_genome = get_human_reference_genome()
    all_transcripts = []
    for name, transcript_exons in tqdm(splitted_transcripts_info):
        gene_name, transcript_name, chromosome, strand = name
        exon_sequences = []
        transcript_sequence = Seq("")
        exons_to_examin = [
            int(k) for k in transcript_exons["Exon Rank in Transcript"].values
        ]
        if strand == "1":
            for exon_number in range(min(exons_to_examin), max(exons_to_examin) + 1):
                exon = transcript_exons[
                    transcript_exons["Exon Rank in Transcript"] == str(exon_number)
                ]
                end = int(exon["3' UTR End"])
                start = int(exon["3' UTR Start"])
                exon_seq = human_genome[chromosome].seq[(start - 1) : end]
                transcript_sequence = transcript_sequence + exon_seq
                exon_sequences.append(exon_seq)
        else:
            for exon_number in reversed(
                range(min(exons_to_examin), max(exons_to_examin) + 1)
            ):
                exon = transcript_exons[
                    transcript_exons["Exon Rank in Transcript"] == str(exon_number)
                ]
                end = int(exon["3' UTR End"])
                start = int(exon["3' UTR Start"])
                exon_seq = human_genome[chromosome].seq[(start - 1) : end]
                transcript_sequence = transcript_sequence + exon_seq
                exon_sequences.append(exon_seq)

        if strand == "-1":
            transcript_sequence = transcript_sequence.reverse_complement()

        transcript_sequence = transcript_sequence.transcribe()
        transcript_record = SeqRecord(
            transcript_sequence,
            id=transcript_name,
            name=gene_name,
            description="3'UTR//{}//{}".format(gene_name, transcript_name),
        )
        all_transcripts.append(transcript_record)

    with open(output_path, "w") as output_handle:
        SeqIO.write(all_transcripts, output_handle, "fasta")


def main():
    # We first want to get a fasta file with every miRNA we want to take into account
    expressed_miRNAs = get_list_of_expressed_miRNAs()

    # Now we want to create a fasta file containing each isomir
    print("get miRNA sequences")
    isomir_fasta_file = "../temp_data/isomir.fasta"
    create_fasta_file_for_isomir(isomir_fasta_file)

    # We also want to create the fasta file for every transcript we are interested in
    print("get 3'UTR sequences")
    transcript_fasta_file = "../temp_data/temp_transcript.fasta"
    get_RNA_of_each_transcript_in_fasta(transcript_fasta_file)

    # Launch miranda
    print("Running miRanda")
    output_miranda = "../results/isomiR_binding_sites_on_protein_coding_3primeUTR.txt"
    run(
        [
            "../miRanda_install/bin/miranda",
            isomir_fasta_file,
            transcript_fasta_file,
            "-out",
            output_miranda,
            "-strict",
            "-quiet",
        ]
    )


########
##Main##
########
if __name__ == "__main__":
    print("Hello world !")
    main()
    print("python script ended succesfully")
