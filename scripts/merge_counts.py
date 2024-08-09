import pandas as pd
import os
import shutil
from concurrent.futures import ThreadPoolExecutor

def read_file(file):
        return pd.read_csv(file, sep="\t", header=0, comment="#",usecols=[6])


#Fusionne les fichiers crée par futureCOunts en une seuls fichier merged_file
def merging_files(summary_file, merged_file, threads=1):

    #Lire de la Datafram à partir du ficher summary
    df_summary = pd.read_csv(summary_file,sep="\t", header=0)

    #Si le fichier n'existe pas le crée en etant la copie du premier fichier
    if not os.path.exists(merged_file) or os.path.getsize(merged_file) == 0:
        shutil.copy(df_summary.iloc[0]["path_counts"],merged_file)

    df_merged = pd.read_csv(merged_file, sep="\t",header=0,comment="#")

    #Extraction les ids et chemins du Dataframe
    paths_bam = df_summary["path_bam"].values.tolist()
    paths_counts = df_summary["path_counts"].values.tolist()

    #Identifier les nouveau fichiers à traiter
    new_files = []
    for i,file in enumerate(paths_bam):
        if file not in df_merged.columns:
            new_files.append(paths_counts[i])

    #Lire les nouveaux fichiers en prallèle
    if new_files:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            data_list = list(executor.map(read_file, new_files))

        #concaténation des données
        new_data = pd.concat(data_list, axis=1, join="outer")
        df_merged = pd.concat([df_merged, new_data],axis=1, join="outer")

        #Ecriture des données fusionnées.
        df_merged.to_csv(merged_file, sep="\t",index=False)