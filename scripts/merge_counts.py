import pandas as pd
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor

def read_file(file, name_col):
    df = pd.read_csv(file, sep="\t", header=0, comment="#",usecols=[6])
    df.columns = [name_col]
    return df

#Fusionne les fichiers crée par futureCOunts en une seuls fichier merged_file
def merging_files(summary_file, merged_file, threads=1):

    #Lire de la Datafram à partir du ficher summary
    df_summary = pd.read_csv(summary_file,sep="\t", header=0)

    #Si le fichier n'existe pas le crée en etant la copie du premier fichier
    if not os.path.exists(merged_file) or os.path.getsize(merged_file) == 0:
        old_id = df_summary.iloc[0]["path_bam"]
        new_id = df_summary.iloc[0]["id"]
        shutil.copy(df_summary.iloc[0]["path_counts"],merged_file)
        subprocess.run(['sed', '-i', f's~{old_id}~{new_id}~g', merged_file])

    df_merged = pd.read_csv(merged_file, sep="\t",header=0,comment="#")

    #Extraction les ids et chemins du Dataframe
    paths_bam = df_summary["path_bam"].values.tolist()
    paths_bam = paths_bam[1:]
    paths_counts = df_summary["path_counts"].values.tolist()
    paths_counts = paths_counts[1:]
    ids = df_summary["id"].values.tolist()
    ids = ids[1:]

    #Identifier les nouveau fichiers à traiter
    new_files = []
    new_ids = []
    for i,file in enumerate(paths_bam):
        if file not in df_merged.columns:
            new_files.append(paths_counts[i])
            new_ids.append(ids[i])

    #Lire les nouveaux fichiers en prallèle
    if new_files:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            data_list = list(executor.map(read_file, new_files, new_ids))

        #concaténation des données
        new_data = pd.concat(data_list, axis=1, join="outer")
        df_merged = pd.concat([df_merged, new_data],axis=1, join="outer")

        #Ecriture des données fusionnées.
        df_merged.to_csv(merged_file, sep="\t",index=False)