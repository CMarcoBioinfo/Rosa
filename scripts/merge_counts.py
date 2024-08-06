import pandas as pd
from concurrent.futures import ThreadPoolExecutor

def read_file(file):
    return pd.read_csv(file, sep="\t", index_col=0, comment="#")


#Fusionne les fichiers crée par futureCOunts en une seuls fichier merged_file
def merging_files(summary_file, merged_file, threads=1):

    #Lire de la Datafram à partir du ficher summary
    df = pd.read_csv(summary_file,sep="\t", index_col=0)

    #Lire le fichier merged existant
    try:
        merged_data = pd.read_csv(merged_file, sep="\t",index_col=0)
    except FileExistsError:
        merged_data = pd.DataFrame()

    #Extraction les ids et chemins du Dataframe
    file_paths = df["path"].tolist()
    ids = df["id"].tolist()

    #Identifier les nouveau fichiers à traiter
    new_files = []
    new_ids = []
    for file, id in zip(file_paths, ids):
        if id not in merged_data.columns:
            new_files.append(file)
            new_ids.append(id)
    
    #Lire les nouveaux fichiers en prallèle
    with ThreadPoolExecutor(max_workers=threads) as executor:
        data_list = list(executor.map(read_file, new_files))

    for i, df in enumerate(data_list):
        df.columns = [new_ids[i]]

    #concaténation des données
    new_data = pd.concat(data_list, axis=1, join="outer")
    merged_data = pd.concat([merged_data, new_data],axis=1, join="outer")

    #Ecriture des données fusionnées.
    merged_data.to_csv(merged_file, sep="\t")