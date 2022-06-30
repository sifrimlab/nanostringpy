import re
import json
import os
import pandas as pd
import itertools
from glob import glob


def getGeneDict(file):
    with open(file, 'r') as meta_file:
        data = json.load(meta_file)
    genes = { target["Probes"][0]["RTS_ID"]:target["DisplayName"] for target in data["Targets"] }
    return genes

def getCountDict(dcc_file, pkc_file, gene_dict = {}):
    if not gene_dict:
        gene_dict = getGeneDict(pkc_file)
    with open(dcc_file, 'r') as count_file:
        gene_counts_dict = {}
        fileContent = count_file.read()
        # r'\[(.*?)\]'
        # print( re.findall(r'\<Code_Summary>(.*?)\<\/Code_Summary\>', fileContent) )
        counts = re.search(r'(?<=<Code_Summary>\n)((.|\n)+)(?=<\/Code_Summary>)', fileContent).group().split("\n")
        for count in counts:
            try:
                rts_id, count = count.split(",")
                gene = gene_dict[rts_id]
                gene_counts_dict[gene] = gene_counts_dict.get(gene, 0) + int( count )
            except:
                continue
    return gene_counts_dict

def createMetadataDf(txt_file):
    df = pd.read_csv(txt_file, sep="\n", header = None)
    r = []
    for i in range(len(df.index)):
        row = df.iloc[i][0]
        if (row.startswith("Sample_ID")):
            h = row.split("\t")
        if (row.startswith("DSP")):
            # print(df.iloc[i+1][0].split("\t"))
            new_row = []
            for i,item in enumerate(row.split("\t")):
                if i == 1:
                    subitems = item.split(" ")
                    for subitem in subitems:
                        new_row.append(subitem)
                else:
                    new_row.append(item)
            r.append(new_row)
    t = pd.DataFrame(data=r, columns=h)
    return t


def getPanelFromSampleID(sample_id, file = "", metadata_df = None):
    if file:
        if metadata_df is not None:
            raise ValueError("Either file or metadata_df should be given, not both.")
        else:
            metadata_df = createMetadataDf(file)
    elif metadata_df is None:
        raise ValueError("Either file or metadata_df should be given.")

    return metadata_df[metadata_df["Sample_ID"] == sample_id]["panel"].iloc[0]
    
def getSampleNameFromSampleID(sample_id, file = "", metadata_df = None):
    if file:
        if metadata_df is not None:
            raise ValueError("Either file or metadata_df should be given, not both.")
        else:
            metadata_df = createMetadataDf(file) if file.endswith(".txt") else pd.read_csv(file)
    elif metadata_df is None:
        raise ValueError("Either file or metadata_df should be given.")

    sample_row = metadata_df[metadata_df["Sample_ID"] == sample_id].iloc[0]

    return f"{sample_row['scan name'].replace(' ', '-')}_{sample_row['segment'].replace(' ', '-')}_{sample_row['group'].replace(' ', '-')}" if sample_row['slide name'] != "No Template Control" else "No Template Control"


def createCountMatrix(dcc_glob_pattern, pkc_file, meta_file):
    gene_dict = getGeneDict(pkc_file)
    metadata_df = createMetadataDf(meta_file) if meta_file.endswith(".txt") else pd.read_csv(meta_file)

    rows_list = []
    header = ["TargetName"]
    for file in glob(dcc_glob_pattern):
        gene_counts_dict = getCountDict(file, pkc_file, gene_dict = gene_dict)
        gene_counts_dict["TargetName"] = getSampleNameFromSampleID(os.path.splitext( os.path.basename(file))[0], metadata_df = metadata_df )
        rows_list.append(gene_counts_dict)

    result_df = pd.DataFrame(rows_list)
    result_df.index = result_df["TargetName"]
    result_df= result_df.drop('TargetName', 1)
    result_df = result_df.transpose()
    result_df.to_csv("count_matrix.csv")

if __name__ == '__main__':
    createCountMatrix("../Raw data files/WTA/DCC/*.dcc", "../Raw data files/WTA/Mm_R_NGS_WTA_v1.0/Mm_R_NGS_WTA_v1.0.pkc", "../Raw data files/WTA/22Apr2021_MsWTA_20210427T2207_LabWorksheet DP.csv")


