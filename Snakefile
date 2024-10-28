#Import dependences python
import subprocess
import os
import math
import shutil
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
    base = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/"
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
        os.makedirs(directory,exist_ok=True)

#Deletes a directory and its contents if is exists
def delete_directory(directory):
    shutil.rmtree(directory, ignore_errors=True)

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
    metadata = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/2-Counts/sample_names.txt"
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
    metadata = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/2-Counts/sample_names.txt"
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
        wd = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/reads1/"
    elif ( read_suffix == "2"):
        wd = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/reads2/"
    for file_ext in [".fastq",".fq"]:
        for file_r in ["_" + read_suffix, "_r" + read_suffix, "_R" + read_suffix, "." + read_suffix, ".r" + read_suffix, ".R" + read_suffix]:
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

def is_gz(file):
    if file.endswith(".gz"):
        return True
    else:
        return False




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

def csv_to_dict(csv):
    df = pd.read_csv(csv,sep = "\t")
    dictionnary = df.groupby("id")["path"].apply(list).to_dict()
    return dictionnary

#Returns a dico of samples
#Scans the directories in the samples directory and lists the 'fastq', '.fq' and 'bam' files.
#Them extracts the identifiers of these files
def dict_samples_directory():
    base = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/"
    reads_directories = ["reads1/","reads2/"]
    format_directories = ["bam/", "reads/"]
    samples = {}

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
                    path = base + directory + sample
                    if ( id not in samples ):
                        samples[id] = [path]
                    else :
                        samples[id].append(path)
    
    #Browse formatting directories 
    for directory in format_directories:
        wd = os.path.abspath(base + directory)
        sample_list = os.listdir(os.path.abspath(wd))

        #Browse file in directories
        for sample in sample_list:
            if ( ".fastq" in sample or ".fq" in sample or ".bam" in sample):
                id = get_id(sample)
                path = base + directory + sample
                samples[id] = path
    
    #Uniq ID

    return samples

samples_directory_date = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/.samples.date"
mergeFile = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/2-Counts/merge.counts"

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



#Obtien le chemin d'entrée des dictionnaire fastq.
def get_inputs(wildcards, mode):
    id = wildcards.id.rsplit("/")[-1]
    if ( mode == 4 ):
        return bam[id]
    else:
        if ( mode == 0 ):
            if ("reads1" in wildcards.id):
                return str(fastq2[id][0])
            elif ("reads2" in wildcards.id):
                return str(fastq2[id][1])
            elif ("reads" in wildcards.id):
                return str(fastq[id])
        if ( mode == 1):
            file = str(fastq2[id][0])
        elif ( mode == 2):
            file = str(fastq2[id][1])
        elif ( mode == 3 ):
            file = str(fastq[id])
        if file.endswith(".gz"):
            return file 
        elif os.path.exists(file):
            return file
        else :
            return file + ".gz"

#Temporaries

tmp_directory = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/.tmp/"
tmp_samples = tmp_directory + "samples/"
tmp_bam = tmp_directory + "bam/"
tmp_reads = tmp_samples + "reads/"
tmp_reads1 = tmp_samples + "reads1/"
tmp_reads2 = tmp_samples + "reads2/"
create_directory_if_not_exists(tmp_directory)
create_directory_if_not_exists(tmp_samples)
create_directory_if_not_exists(tmp_reads)
create_directory_if_not_exists(tmp_reads1)
create_directory_if_not_exists(tmp_reads2)
create_directory_if_not_exists(tmp_bam)

#Remplit les différents dictionnaire fastq, fastq2 et bam et renvoie le liste des id des que le snakemake se lance normalement.
def all_samples(inputs_dict, fastq2, fastq, bam):
    for key, item in inputs_dict.items():
        if (len(item) == 2):
            if ( (os.path.exists(item[0]) and os.path.exists(item[1])) or (os.path.exists(item[0] + ".gz") and os.path.exists(item[1] + ".gz")) or (os.path.exists(item[0]) and os.path.exists(item[1] + ".gz")) or (os.path.exists(item[0] + ".gz") and os.path.exists(item[1])) ):
                if ( ((".fastq" in item[0] or ".fq" in item[0])) and ((".fastq" in item[1] or ".fq" in item[1])) ) :
                    read1 = str(tmp_reads1) + str(key) + ".pre"
                    read2 = str(tmp_reads2) + str(key) + ".pre"
                    fastq2[str(key)] = [item[0],item[1]]
                    os.system(f"touch {read1} {read2}")

        elif (len(item)== 1 and (os.path.exists(item[0]) or os.path.exists(item[0] + "gz") )):
            if( ".bam" in item[0]):
                bam_file = str(tmp_bam) + str(key) + ".pre"
                bam[str(key)] = item[0]
                os.system(f"touch {bam_file}")
            elif ( ".fastq" in item[0] or ".fq" in item[0]):
                read = str(tmp_reads) + str(key) + ".pre"
                fastq[str(key)] = item[0]
                os.system(f"touch {read}")
 

    all_samples = set(list(fastq.keys()) + list(fastq2.keys()) + list(bam.keys()))
    return all_samples


#Variables

fastq2 = {}
fastq = {}
bam = {}


csv = config["GENERAL"]["DATA_INPUTS"]["SAMPLES_FILE"]
if(check_value(csv)) :

    dict_csv = csv_to_dict(csv)
    all_samples = all_samples(dict_csv, fastq2, fastq, bam)


else:
    dict_directory = dict_samples_directory()
    all_samples = all_samples(dict_directory, fastq2, fastq, bam)


#     for key, item in dictionnary.items():
#         if (len(item) == 2):
#             if ( (os.path.exists(item[0]) and os.path.exists(item[1])) or (os.path.exists(item[0] + ".gz") and os.path.exists(item[1] + ".gz")) ):
#                 if ( ((".fastq" in item[0] or ".fq" in item[0])) and ((".fastq" in item[1] or ".fq" in item[1])) ) :
#                     read1 = tmp_reads1 + key + ".pre"
#                     read2 = tmp_reads2 + key + ".pre"
#                     fastq2[key] = [item[0],item[1]]
#                     os.system(f"touch {read1} {read2}")


#         elif (len(item)== 1 and (os.path.exists(item[0]) or os.path.exists(item[0] + "gz") )):
#             if( ".bam" in item[0]):
#                 bam_file = tmp_bam + key + ".pre"
#                 bam[key] = item[0]
#                 os.system(f"touch {bam_file}")
#             elif ( ".fastq" in item[0] or ".fq" in item[0]):
#                 read = tmp_read + key + ".pre"
#                 fastq[key] = item[0]
#                 os.system(f"touch {read}")

#     all_samples = set(list(fastq.keys()) + list(fastq2.keys()) + list(bam.keys()))

#     

# else :
#     all_samples = list_samples()

#     include: "modules/formatting_directories.smk"


summary_file = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/2-Counts/sample_names.txt"
if(os.path.exists(summary_file)) :
    df_summary = pd.read_csv(summary_file,sep="\t", header=0)
    ids = df_summary["id"].values.tolist()


#Obtenir l'heure actuelle
current_time = time.localtime()

#Formater la date et l'heure comme identifiant unique
unique_id = time.strftime("%Y%m%d%H%M%S", current_time)
 


#Calling Snakemake module
include: "modules/formatting_file.smk"
include: "modules/quality_control_fastq.smk"
include: "modules/spliceLaucher.smk"
#include: "modules/featureCounts.smk"
#include: "modules/hisat2.smk"
include: "modules/quality_control_bam.smk"


# rule all:
#     input:
#         mergeFile = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/2-Counts/merge.counts"

# rule all:
#     input:
#         # read1 =expand(os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples/{reads}_1.fastq.gz"),reads= all_samples),
#         # read2 =expand(os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples/{reads}_2.fastq.gz"),reads= all_samples),
#         #bam = expand(config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam.bai",reads= all_samples),
#         directory_data = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/multiqc/" + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] +"_fastq_raw_" + unique_id + "_data/",
#         html =config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/multiqc/" + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"]  + "_fastq_raw_" + unique_id + ".html",
#         directory_data2 = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/multiqc/" + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] +"_fastq_trimmed_" + unique_id + "_data/",
#         html2 = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"]  + "/2-processed_data/quality_control/multiqc/" + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"]  + "_fastq_trimmed_" + unique_id + ".html",
#         directory_data_bam = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/multiqc/" + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] +"_bam_" + unique_id + "_data/",
#         html_bam = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/multiqc/" + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"]  + "_bam_" + unique_id + ".html",


rule all:
    input:
        reference = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"],
        SA = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/SA",
        SAindex = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/SAindex"
rule delete_tmp:
    input: 
        rules.all.output

    params:
        tmp_directory = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/.tmp/"
    
    run:
        delete_directory(tmp_directory)