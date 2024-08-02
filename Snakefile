#Import dependences python
import subprocess
import os
import sys


#Function use

#Create directory
def create_directory_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory,exsit_ok=True)

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
    reads_directories = ["reads1","reads2"]
    format_directories = ["bam/", "reads"]
    samples = {"id": []}

    #Browse reads directories
    for directory in reads_directories:
        wd = base + directory
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
        wd = base + directory
        sample_list = os.listdir(os.path.abspath(wd))

        #Browse file in directories
        for sample in sample_list:
            if ( ".fastq" in sample or ".fq" in sample or ".bam" in sample):
                id = get_id(sample)
                samples["id"].append(id)
    
    #Uniq ID
    all_samples = set(samples["id"])

    return all_samples






#Variable

all_samples = list_samples()


#Calling Snakemake module

include: "modules/formatting.smk"
include: "modules/hisat2.smk"
include: "modules/featureCounts.smk"


rule all:
    input:
        controlFiles = expand(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/.Control/{all_samples}.txt",all_samples=all_samples)