import os
import pandas as pd
from filelock import FileLock
from concurrent.futures import ThreadPoolExecutor

def read_csv_file(csv, all_columns=False):
    if all_columns:
            return pd.read_csv(csv, sep="\t", index_col=0, comment="#")
    return pd.read_csv(csv, sep="\t", index_col=0, usecols=[0, 6], comment="#")

def merging_files(merged_file, count_file, control_file, threads=1):
    
    #Create Filelock object
    directory = os.path.dirname(merged_file)
    file_name = os.path.basename(merged_file)
    lock_file_name = "." + file_name + ".lock"
    lock_file = os.path.join(directory,lock_file_name)
    lock = FileLock(lock_file)

    with lock:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            #Read count file
            all_columns = not os.path.exists(merged_file)
            executor_count_file = executor.submit(read_csv_file, count_file, all_columns)
            count_file_df = executor_count_file.result()

            #Check if merged file exists
            if os.path.exists(merged_file) and os.path.getsize(merged_file) > 0 :
                #Read merged file
                executor_merged_file = executor.submit(pd.read_csv, merged_file, sep="\t", index_col=0)
                merged_file_df = executor_merged_file.result()
            
                #Add new count file
                merged_df = pd.concat([merged_file_df,count_file_df], axis=1, join="inner")

            else:
                #Use count file as initial merged file
                merged_df = count_file_df

        #Save merge file
        merged_df.to_csv(merged_file, sep="\t")

        with open(control_file, "a") as file:
            file.write(f"{count_file}\n")