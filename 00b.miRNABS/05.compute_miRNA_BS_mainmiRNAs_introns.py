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


def create_fasta_file_for_miRNA(path_fasta_file, miRNAs_names):
    # First wa want to be abble to find the sequence of each miRNAs, so we open mirbase20 in a fasta format and filter it
    sequence_for_each_miRNA = {}

    for record in SeqIO.parse(
        "/pasteur/projets/policy01/evo_immuno_pop/ERCPilot_SharedDBs/mirbase20/mature.fa",
        "fasta",
    ):
        if record.id in miRNAs_names:
            if record.id in sequence_for_each_miRNA.keys():
                if sequence_for_each_miRNA[record.id] != record.seq:
                    miRNAs_names.remove(record.id)
                    print("removed {}".format(record.id))
            else:
                sequence_for_each_miRNA[record.id] = record.seq

    # conversion en list
    list_of_record_to_output = [
        SeqRecord(s, id=n, description="", name="")
        for n, s in sequence_for_each_miRNA.items()
    ]

    with open(path_fasta_file, "w") as output_handle:
        SeqIO.write(list_of_record_to_output, output_handle, "fasta")


def get_human_reference_genome():
    fasta_human_genome = "/pasteur/projets/policy01/evo_immuno_pop/ERCPilot_SharedDBs/HumanGenomeGRCh37/GRCh37_dna_sm_ref.fa"
    record_human = SeqIO.to_dict(SeqIO.parse(fasta_human_genome, "fasta"))
    return record_human


def setdiff_range(list_first, list_second):
    start_intron =  int(list_first[0])
    end_intron =  int(list_first[1])
    start_exon =  int(list_second[0])
    end_exon =  int(list_second[1])

    #First e treat the case were there is no intersection
    if (end_exon < start_intron) or (start_exon > end_intron) :
        #No intersection
        return([[start_intron, end_intron]])

    #We treat the case where the introns completely inside the exon
    if (start_exon <= start_intron) and (end_exon >= end_intron) :
        return([])

    #we treat the case where the exon in inside the intron
    if (start_exon > start_intron) and (end_exon < end_intron) :
        return([[start_intron, start_exon-1], [end_exon +1, end_intron]])

    #the exon intersect partially with the left of the intron
    if (start_exon <= start_intron) and (end_exon< end_intron) :
        return([[end_exon+1, end_intron]])

    if (start_exon>start_intron) and (end_exon>= end_intron):
        return([[start_intron, start_exon -1 ]])

    print("PROBLEM ! CASE NOT CATCHED")

def get_RNA_of_introns(output_path):
    print("getting introns coordinates")
    exons_informations = pd.read_csv(
        "/pasteur/projets/policy01/evo_immuno_pop/Maxime/Evo_Immuno_pop_data/ExonCoordinates_hg37_ens70.txt",
        dtype="str",
        sep="\t",
    )
    exons_informations.sort_values(by=["Chromosome Name"], inplace=True)
    # We only want to keep lines where chromosome is in the "chromosome_to_keep" list
    chromosome_to_keep = [str(k) for k in range(1, 23)]
    chromosome_to_keep.append("X")
    exons_informations = exons_informations[
        exons_informations["Chromosome Name"].isin(chromosome_to_keep)
    ]

    exons_informations = exons_informations[["Ensembl Gene ID", "Ensembl Transcript ID", "Chromosome Name", "Exon Chr Start (bp)", "Exon Chr End (bp)", "Strand"]]

    exons_informations.rename(columns={"Ensembl Gene ID" : "Ensembl_Gene_ID",
                                        "Ensembl Transcript ID" : "Ensembl_Transcript_ID",
                                        "Chromosome Name" : "Chromosome_Name",
                                        "Exon Chr Start (bp)" : "Exon_Chr_Start",
                                        "Exon Chr End (bp)" : "Exon_Chr_End"}, inplace=True)
    all_exons_per_gene = exons_informations.groupby(["Ensembl_Gene_ID", "Chromosome_Name", "Strand"])
    all_introns = []
    for name, gene_exons in tqdm(all_exons_per_gene) :
        gene_name, chromosome, strand = name
        introns = [ [min( [ int(k) for k in gene_exons["Exon_Chr_Start"].values]),
                     max(  [int(k) for k in gene_exons["Exon_Chr_End"].values])] ]
        for gene_exon in gene_exons.itertuples() :
            exon = [int(gene_exon.Exon_Chr_Start), int(gene_exon.Exon_Chr_End)]
            new_introns = []
            for intron in introns :
                new_intron_s = setdiff_range(intron, exon)
                for k in new_intron_s :
                    new_introns.append(k)
            introns = new_introns.copy()
        intron_start = [k[0] for k in introns]
        intron_end = [k[1] for k in introns]
        introns_df = pd.DataFrame({'gene_name' : gene_name,
                        'chromosome' : chromosome,
                        'strand' : strand,
                       'intron_number' : range(1, len(intron_end) +1 ),
                       'intron_start' : intron_start,
                        'intron_end' : intron_end})
        all_introns.append(introns_df)
    all_introns_df = pd.concat(all_introns)

    #Now that we have coordinate (and strand), we create the fasta file
    all_seq_record = []
    human_genome = get_human_reference_genome()
    print("getting introns sequence")
    for intron in tqdm(all_introns_df.itertuples()) :
        if intron.strand == "+" :
            intron_seq = human_genome[intron.chromosome].seq[int(intron.intron_start) -1 : int(intron.intron_end)]
            intron_seq = intron_seq.transcribe()
        else :
            intron_seq = human_genome[intron.chromosome].seq[int(intron.intron_start) -1 : int(intron.intron_end)].reverse_complement()
            intron_seq = intron_seq.transcribe()
        intron_record = SeqRecord(
            intron_seq,
            id="{}_intron_{}".format(intron.gene_name, intron.intron_number),
            name=gene_name,
            description="intron_{}_{}_{}_{}".format(intron.gene_name, intron.chromosome, int(intron.intron_start), int(intron.intron_end)),
        )
        all_seq_record.append(intron_record)
    with open(output_path, "w") as output_handle:
        SeqIO.write(all_seq_record, output_handle, "fasta")
def main():
    # We first want to get a fasta file with every miRNA we want to take into account
    expressed_miRNAs = get_list_of_expressed_miRNAs()
    #
    # # Now we want to create a fasta file containing each miRNA
    print("get miRNA sequences")
    miRNA_fasta_file = "../temp_data/miRNA2.fasta"
    create_fasta_file_for_miRNA(miRNA_fasta_file, expressed_miRNAs)

    # We also want to create the fasta file for every transcript we are interested in
    # print("get introns sequences")
    transcript_fasta_file = "../temp_data/temp_introns.fasta"
    get_RNA_of_introns(transcript_fasta_file)


    # Launch miranda
    print("Running miRanda")
    output_miranda = "../results/miRNA_binding_sites_on_protein_coding_introns.txt"
    run(
        [
            "../miRanda_install/bin/miranda",
            miRNA_fasta_file,
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
