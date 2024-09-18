#Import dependences python
import subprocess
import os
import time
from pathlib import Path
import sys
import pandas as pd
import datetime
import time
from inotify_simple import INotify, flags
import fcntl

#Function use

#Obtient la date de modification d'un path
def get_date(path):
    return datetime.datetime.fromtimestamp(os.path.getmtime(path))

#Ecris la date d'un fichier ou dossier dans un fichier
def write_date(path,date):
    with open(path,"w") as file:
        file.write(date.isoformat())

#Lire le fichier de date
def read_date(path):
    with open(path, "r") as file:
        date_str = file.read().strip()
        return datetime.datetime.fromisoformat(date_str)

#Renvoie la date la plus courrente des dossiers.

def directories_recent_data():
    base = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/"
    reads_directories = ["reads1/","reads2/","bam","reads"]
    date = None
    for directory in reads_directories:
        wd = os.path.abspath(base + directory)
        tmp_date = get_date(wd)
        if date is None or tmp_date > date:
            date = tmp_date
    return date


#Create directory
def create_directory_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory,exsit_ok=True)

#Lock un fichier
def lock_file(file_path):
    if not os.path.exists(file_path):
        open(file_path, "w").close()
    lock_file = open(file_path, "r+")
    fcntl.flock(lock_file, fcntl.LOCK_EX)
    return lock_file

#Unlock un fichier
def unclock_file(lock_file):
    fcntl.flock(lock_file, fcntl.LOCK_UN)
    lock_file.close()

#Ajoute un sample au fichier metadata
def add_file(name,path_bam,path_counts):
    metadata = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/sample_names.txt"
    file_lock= lock_file(metadata)
    if not os.path.exists(metadata) or os.path.getsize(metadata) == 0:
        df = pd.DataFrame(columns=["id", "path_bam", "path_counts"])

    else:
        df = pd.read_csv(metadata,sep="\t",header=0,dtype={"id":str, "path_bam":str, "path_counts":str})

    if name not in df["id"].values:
        new_row = pd.DataFrame([{"id" : name, "path_bam": path_bam, "path_counts": path_counts}])
        df = pd.concat([df,new_row], ignore_index=True)
        df.to_csv(metadata, sep="\t", index=False)
    unclock_file(file_lock)

#Return list of available count files
def get_available_files():
    metadata = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/sample_names.txt"
    if os.path.exists(metadata):
        df = pd.read_csv(metadata, sep = "\t")
        files = df["path"].tolist()
    else:
        files = []
    return files




#Check if varaible is None or empty
def check_value (var):
    if var is None:
        return False
    elif type(var) is int:
        return True
    elif type(var) is str:
        if len(var.strip()) == 0:
            return False
        else:
            return True

#Returns the list of reads (sample.fastq.gz)
def reads(wcs, read_suffix, filter_fastq=False):
    name = wcs.reads
    if (read_suffix == "1"):
        wd = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/reads1/"
    elif ( read_suffix == "2"):
        wd = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/reads2/"
    for file_ext in [".fastq",".fq"]:
        for file_r in ["_" + read_suffix, "_r" + read_suffix, "_R" + read_suffix, "." + read_suffix, ".r" + read_suffix, "R" + read_suffix]:
            file = wd + name + file_r +file_ext
            if (os.path.exists(file) or os.path.exists(file + ".gz")):
                read = wd + name + file_r + file_ext + ".gz"
    if filter_fastq:
        read = wd + name + "_" + read_suffix + ".filter.fastq.gz"
    read = os.path.abspath(read)
    return read



#Extract file ID
#If the file name ends with '.gz', it deletes the last two extension (.fq.gz)
#Else, it deletes the last extension
def get_id(file):
    if file.endswith(".gz"):
        id = file.rsplit(".", 2)[0]
    else:
        id = file.rsplit(".", 1)[0]
    return id


#Checks the suffix and return ID without the suffix
def check_suffix(id):
    ext_underscore = id.rsplit("_", 1)
    ext_point = id.rsplit(".", 1)
    suffixes = ["1", "2", "r1", "r2", "R1", "R2"]
    if ( len(ext_underscore) == 2 and ext_underscore[1] in suffixes ):
        return ext_underscore[0]
    if ( len(ext_point) == 2 and ext_point[1] in suffixes ):
        return ext_point[0]
    return None




#Returns a list of samples
#Scans the directories in the samples directory and lists the 'fastq', '.fq' and 'bam' files.
#Them extracts the identifiers of these files
def list_samples():
    base = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/"
    reads_directories = ["reads1/","reads2/"]
    format_directories = ["bam/", "reads/"]
    samples = {"id": []}

    #Browse reads directories
    for directory in reads_directories:
        wd = os.path.abspath(base + directory)
        sample_list = os.listdir(os.path.abspath(wd))
    
        #Browse file in directories
        for sample in sample_list:
            if ( ".fastq" in sample or ".fq" in sample ):
                raw_id = get_id(sample)
                id = check_suffix(raw_id)
                if ( id is not None ):
                    samples["id"].append(id)
    
    #Browse formatting directories 
    for directory in format_directories:
        wd = os.path.abspath(base + directory)
        sample_list = os.listdir(os.path.abspath(wd))

        #Browse file in directories
        for sample in sample_list:
            if ( ".fastq" in sample or ".fq" in sample or ".bam" in sample):
                id = get_id(sample)
                samples["id"].append(id)
    
    #Uniq ID
    all_samples = set(samples["id"])
    return all_samples

samples_directory_date = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/.samples.date"
mergeFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/merge.counts"

if (not os.path.exists(samples_directory_date)):
    initial_directory_date = directories_recent_data()
    write_date(samples_directory_date,initial_directory_date)
else:
    initial_directory_date = read_date(samples_directory_date)

current_directory_date = directories_recent_data()
if (os.path.exists(mergeFile)):
    if (current_directory_date.timestamp() != initial_directory_date.timestamp()):
        current_time = time.time()
        os.utime(mergeFile,(current_time, current_time))
        write_date(samples_directory_date,current_directory_date)


#Variables
summary_file = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/sample_names.txt"
if(os.path.exists(summary_file)) :
    df_summary = pd.read_csv(summary_file,sep="\t", header=0)
    ids = df_summary["id"].values.tolist()

all_samples = list_samples()


#Obtenir l'heure actuelle
current_time = time.localtime()

#Formater la date et l'heure comme identifiant unique
unique_id = time.strftime("%Y%m%d%H%M%S", current_time)
 


#Calling Snakemake module

include: "modules/formatting.smk"
include: "modules/quality_control_fastq.smk"
include: "modules/hisat2.smk"
include: "modules/featureCounts.smk"


# rule all:
#     input:
#         mergeFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/merge.counts"


rule all:
    input:
        # read1 =expand(os.path.abspath(config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/samples/{reads}_1.fastq.gz"),reads= all_samples),
        # read2 =expand(os.path.abspath(config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/samples/{reads}_2.fastq.gz"),reads= all_samples),
        #bam = expand(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam.bai",reads= all_samples),
        directory_data = directory(config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/multiqc/" + config["PARAMS"]["GENERAL"]["PREFIX"] +"_fastq_raw_" + unique_id + "_data/"),
        html =config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/multiqc/" + config["PARAMS"]["GENERAL"]["PREFIX"]  + "_fastq_raw_" + unique_id + ".html",
        directory_data2 = directory(config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/multiqc/" + config["PARAMS"]["GENERAL"]["PREFIX"] +"_fastq_trimmed_" + unique_id + "_data/"),
        html2 = config["DATA_INPUT"]["WORKING_DIRECTORY"]  + "/2-processed_data/quality_control/multiqc/" + config["PARAMS"]["GENERAL"]["PREFIX"]  + "_fastq_trimmed_" + unique_id + ".html"
