#############################################################################
##We want to take the miRanda output and create a file with the interaction##
##Between miRNA and transcript, without keeping into account the location  ##
##of the binding                                                           ##
#############################################################################

import pandas as pd


def get_data_frame_basis():
    infile = "../results/common_isomiR_binding_sites_on_protein_coding_3primeUTR.txt"
    miRNAs = []
    transcripts = []
    miRNABS_start_pos = []
    miRNABS_length = []
    per1 = []
    per2 = []
    with open(infile, "r") as handle:
        for l in handle:
            if (">" == l[0]) and (">" != l[1]):
                working_line = l.strip().split()
                miRNAs.append(working_line[0][1:])
                transcripts.append(working_line[1])
                miRNABS_start_pos.append(working_line[6])
                miRNABS_length.append(working_line[8])
                per1.append(working_line[9])
                per2.append(working_line[10])

    interaction_dictionnary = {}
    interaction_dictionnary["isomirID"] = miRNAs
    interaction_dictionnary["transcripts"] = transcripts
    interaction_dictionnary["isomiRNABS_start_pos"] = miRNABS_start_pos
    interaction_dictionnary["isomiRNABS_length"] = miRNABS_length
    interaction_dictionnary["per1"] = per1
    interaction_dictionnary["per2"] = per2

    interaction_df = pd.DataFrame(interaction_dictionnary)
    del interaction_dictionnary
    return interaction_df


def add_ensembl_gene_ID_and_gene_name(interact_df):
    print("loading genes data")
    transcripts_info = pd.read_table(
        "/Volumes/evo_immuno_pop/Maxime/Evo_Immuno_pop_data/ExonCoordinates_hg37_ens70.txt",
        dtype="str",
    )
    transcripts_info = transcripts_info[
        ["Ensembl Gene ID", "Ensembl Transcript ID", "Associated Gene Name"]
    ]
    transcripts_info = transcripts_info.drop_duplicates()
    transcripts_info.index = transcripts_info["Ensembl Transcript ID"]

    interact_df["EnsemblGeneID"] = interact_df["transcripts"].map(
        transcripts_info["Ensembl Gene ID"]
    )
    interact_df["GeneNames"] = interact_df["transcripts"].map(
        transcripts_info["Associated Gene Name"]
    )

    return interact_df


def main():
    print("loading interaction data")
    interaction_df = get_data_frame_basis()

    interaction_df = add_ensembl_gene_ID_and_gene_name(interaction_df)
    interaction_df.to_csv(
        "../results/common_isomiR_binding_sites_on_protein_coding_3primeUTR_simplified.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    print("Hello World!")
    main()
    print("Python script succesfully ended")
