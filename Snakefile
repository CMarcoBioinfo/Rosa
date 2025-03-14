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
import gc
import psutil


total_ram = psutil.virtual_memory().total
#available_ram = psutil.virtual_memory().available
mem_total_mb = int(total_ram / (1024 * 1024) * 0.8)  # On réserve 80% de la mémoire pour le système

def calculate_genome_length(fasta_file):
    genome_length = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome_length += len(line.strip())
    return genome_length

# cores = min(psutil.cpu_count(), workflow.cores)
#Vérification de l'option dag et rulegraph
args = sys.argv
dag_mode = "--dag" in args
rulegraph_mode = "--rulegraph" in args
filegraph_mode = "--filegraph" in args


if "SNAKEMAKE_PRINT" not in os.environ:
    os.environ["SNAKEMAKE_PRINT"] = "false"

def print_once(message):
    if not dag_mode and not rulegraph_mode and not filegraph_mode:
        if os.getenv("SNAKEMAKE_PRINT") == "false":
            print(message)

if "MAX_CORES" not in os.environ:
    bool_cores = False
    if "--cores" in args:
        index_cores = sys.argv.index("--cores")
        os.environ["MAX_CORES"] = sys.argv[index_cores]
        bool_cores = True

    if "-c" in args :
        index_cores = sys.argv.index("-c") + 1
        os.environ["MAX_CORES"] = sys.argv[index_cores]
        bool_cores = True

    if "--jobs" in args :
        index_cores = sys.argv.index("--jobs") + 1
        os.environ["MAX_CORES"] = sys.argv[index_cores]
        bool_cores = True

    if "-j" in args :
        index_cores = sys.argv.index("-j") + 1
        os.environ["MAX_CORES"] = sys.argv[index_cores]
        bool_cores = True

    if not bool_cores:
        os.environ["MAX_CORES"] = "1"
    
max_cores = int(os.environ["MAX_CORES"])
print_once(f"Valeur de max_cores : {max_cores}")

#Obtenir l'heure actuelle
if os.getenv("SNAKEMAKE_PRINT") == "false":
    current_time = time.localtime()
     #Formater la date et l'heure comme identifiant unique
    unique_id = time.strftime("%Y%m%d%H%M%S", current_time)
    os.environ["UNIQUE_ID"] = unique_id

unique_id = os.environ["UNIQUE_ID"]



#Path data
print_once("Vérification des paramètres")

working_directory = config["GENERAL"].get("WORKING_DIRECTORY")
if not working_directory:
    working_directory = "data"
print_once(f"Working directory : {working_directory} ...... OK")

prefix = config["GENERAL"].get("PREFIX")
if not prefix:
    prefix = "run_" + unique_id
print_once(f"Prefix: {prefix} ...... OK")

path_qc = working_directory + "/3-Quality_control/"
path_results = working_directory + "/4-Results/" + prefix


# #Function use
# #Obtient la date de modification d'un path
# def get_date(path):
#     return datetime.datetime.fromtimestamp(os.path.getmtime(path))


# #Ecris la date d'un fichier ou dossier dans un fichier
# def write_date(path,date):
#     with open(path,"w") as file:
#         file.write(date.isoformat())

# #Lire le fichier de date
# def read_date(path):
#     with open(path, "r") as file:
#         date_str = file.read().strip()
#         return datetime.datetime.fromisoformat(date_str)

#Renvoie la date la plus courrente des dossiers.

# def directories_recent_data():
#     base = working_directory + "/1-raw_data/samples/"
#     reads_directories = ["reads1/","reads2/","bam","reads"]
#     date = None
#     for directory in reads_directories:
#         wd = os.path.abspath(base + directory)
#         tmp_date = get_date(wd)
#         if date is None or tmp_date > date:
#             date = tmp_date
#     return date


def memory_release():
    gc.collect()

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
    metadata = working_directory + prefix + "/2-Counts/sample_names.txt"
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
    metadata = working_directory + prefix + "/2-Counts/sample_names.txt"
    if os.path.exists(metadata):
        df = pd.read_csv(metadata, sep = "\t")
        files = df["path"].tolist()
    else:
        files = []
    return files


#Check if variable is None or empty
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
        wd = working_directory + "/1-raw_data/samples/reads1/"
    elif ( read_suffix == "2"):
        wd = working_directory + "/1-raw_data/samples/reads2/"
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

#Vérifie si l'éxcutable existe
def check_excutable(executable): 
    if shutil.which(executable) is None:
        print_once(f"Erreur : L'exécutable {executable} est manquant.")
        sys.exit(1)
    else : 
        print_once(f"{executable} ...... OK")
    

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


#Returns a dico of samples
#Scans the directories in the samples directory and lists the 'fastq', '.fq' and 'bam' files.
#Them extracts the identifiers of these files
def dict_samples_directory():
    base = working_directory + "/1-raw_data/samples/"
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

# samples_directory_date = working_directory + "/1-raw_data/samples/.samples.date"
# mergeFile = working_directory + prefix + "/2-Counts/merge.counts"

# if (not os.path.exists(samples_directory_date)):
#     initial_directory_date = directories_recent_data()
#     write_date(samples_directory_date,initial_directory_date)
# else:
#     initial_directory_date = read_date(samples_directory_date)

# current_directory_date = directories_recent_data()
# if (os.path.exists(mergeFile)):
#     if (current_directory_date.timestamp() != initial_directory_date.timestamp()):
#         current_time = time.time()
#         os.utime(mergeFile,(current_time, current_time))
#         write_date(samples_directory_date,current_directory_date)



#Obtien le chemin d'entrée des dictionnaire fastq.
def get_inputs(wildcards, mode):
    id = wildcards.id.rsplit("/")[-1]
    if ( mode == 4 ):
        return bam[id]
    else:
        if ( mode == 0 ):
            if ("reads1" in wildcards.id):
                file = str(fastq2[id][0])
            elif ("reads2" in wildcards.id):
                file = str(fastq2[id][1])
            elif ("reads" in wildcards.id):
                file = str(fastq[id])
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

tmp_directory = working_directory + "/.tmp/"
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
# def all_samples(inputs_dict, fastq2, fastq, bam, length_fastp = None):
#     for key, item in inputs_dict.items():
#         if (len(item) == 2):
#             if ( (os.path.exists(item[0]) and os.path.exists(item[1])) or (os.path.exists(item[0] + ".gz") and os.path.exists(item[1] + ".gz")) or (os.path.exists(item[0]) and os.path.exists(item[1] + ".gz")) or (os.path.exists(item[0] + ".gz") and os.path.exists(item[1])) ):
#                 if ( ((".fastq" in item[0] or ".fq" in item[0])) and ((".fastq" in item[1] or ".fq" in item[1])) ) :
#                     if (length_fastp) :
#                         name_key = str(key) + "." + str(length_fastp) + "bp"
#                     else :
#                         name_key = str(key)
#                     read1 = str(tmp_reads1) + name_key + ".pre"
#                     read2 = str(tmp_reads2) + name_key + ".pre"
#                     fastq2[name_key] = [item[0],item[1]]
#                     if (not os.path.exists(read1)):
#                         os.system(f"touch {read1}")
#                     if (not os.path.exists(read2)):
#                         os.system(f"touch {read2}")

#         elif (len(item)== 1 and (os.path.exists(item[0]) or os.path.exists(item[0] + "gz") )):
#             if( ".bam" in item[0]):
#                 if (length_fastp) :
#                     name_key = str(key) + "." + str(length_fastp) + "bp"
#                 else :
#                     name_key = str(key)
#                 bam_file = str(tmp_bam) + name_key + ".pre"
#                 bam[name_key] = item[0]
#                 if (not os.path.exists(bam_file)):
#                     os.system(f"touch {bam_file}")
#             elif ( ".fastq" in item[0] or ".fq" in item[0]):
#                 if (length_fastp) :
#                     name_key = str(key) + "." + str(length_fastp) + "bp"
#                 else :
#                     name_key = str(key)
#                 read = str(tmp_reads) + name_key + ".pre"
#                 fastq[name_key] = item[0]
#                 if (not os.path.exists(read)):
#                     os.system(f"touch {read}")
 

#     all_samples = set(list(fastq.keys()) + list(fastq2.keys()) + list(bam.keys()))
#     return all_samples



def process_sample(dict_csv, length_fastp=None):
    all_samples = set()
    for index, row in dict_csv.iterrows():
        id_sample = row['id']
        read1 = row['path_read1']
        read2 = row['path_read2']
        if (os.path.exists(read1) and os.path.exists(read2)) or \
           (os.path.exists(read1 + ".gz") and os.path.exists(read2 + ".gz")) or \
           (os.path.exists(read1) and os.path.exists(read2 + ".gz")) or \
           (os.path.exists(read1 + ".gz") and os.path.exists(read2)):
            if (".fastq" in read1 or ".fq" in read1) and (".fastq" in read2 or ".fq" in read2):
                name_key = f"{id_sample}.{length_fastp}bp" if length_fastp else f"{id_sample}"
                read1_path = os.path.join(tmp_reads1, name_key + ".pre")
                read2_path = os.path.join(tmp_reads2, name_key + ".pre")
                if not os.path.exists(read1_path):
                    os.system(f"touch {read1_path}")
                if not os.path.exists(read2_path):
                    os.system(f"touch {read2_path}")
                all_samples.add(name_key)

    return all_samples


use_trimming = config["USAGE"].get("TRIMMING")
if use_trimming:
    length_fastp = config["TRIMMING"].get("LENGTH")
    if not length_fastp:
        length_fastp = int(100)
else :
    length_fastp = None

csv = config["GENERAL"].get("SAMPLES_FILE")
all_samples = {}
if(check_value(csv)) :
    dict_csv = pd.read_csv(csv,sep = "\t")
    all_samples = process_sample(dict_csv, length_fastp)

#Vérification des parametre config.yaml

list_inputs = []

#Vérification des executables et paramètres généraux.
if not all_samples:
    print_once(f"Erreur : Aucun patient fournit")
    sys.exit(1)
else:
    print_once("Patients ...... OK")

genome = config["GENERAL"].get("GENOME")
if not os.path.isfile(genome):
    genome = working_directory + "/1-raw_data/reference/" + genome

if not os.path.isfile(genome):
        print_once(f"Erreur : Le fichier génome : {genome} n'existe pas.")
        sys.exit(1)
else : 
    print_once(f"Fichier {genome} ...... OK")
    genome = os.path.abspath(genome)
    name_genome = genome.rsplit(".",1)[0].rsplit("/",1)[1]


pigz = config["DEPENDANCES"]["FORMATING"].get("PIGZ")
if not pigz:
    pigz = "pigz"
check_excutable(pigz)

include: "modules/compress_fastq.smk"
print_once("Module compress_fastq ...... OK")

#Vérification des executable de trimming.
if use_trimming:
    fastp = config["DEPENDANCES"]["FORMATING"].get("FASTP")
    if not fastp:
        fastp = "fastp"
    check_excutable(fastp)
    path_fastq = working_directory + "/2-processed_data/trimmed/fastq/"
    path_bam = working_directory + "/2-processed_data/trimmed/BAM/"
    
    include: "modules/trimming_fastq.smk"
    print_once("Module trimming_fastq ...... OK")

 
else:
    path_fastq = working_directory + "/1-raw_data/fastq/"
    path_bam = working_directory + "/2-processed_data/not_trimmed/BAM/"



#Vérification des executables de qualité controle
use_qc = config["USAGE"].get("QC")
if use_qc:
    fastqc = config["DEPENDANCES"]["QC"].get("FASTQC")
    if not fastqc:
        fastqc = "fastqc"
    check_excutable(fastqc)
    multiqc = config["DEPENDANCES"]["QC"].get("MULTIQC")
    if not multiqc:
        multiqc = "multiqc"
    check_excutable(multiqc)

    include: "modules/quality_control_fastq.smk"
    print_once("Module quality_control_fastq ...... OK")
    #include: "modules/quality_control_bam.smk"
    #print_once("Module quality_control_bam ...... OK")

    #ajout des fichiers de sorties
    directory_data_raw = path_qc + "multiqc/fastq_raw/" + prefix +"_" + unique_id + "_data/"
    html_raw = path_qc + "multiqc/fastq_raw/" + prefix + "_" + unique_id + ".html"
    list_inputs.append(directory_data_raw)
    list_inputs.append(html_raw)
    if use_trimming:
        directory_data_trimmed = path_qc + "multiqc/fastq_trimmed/" + prefix +"_" + unique_id + "_data/"
        html_trimmed = path_qc + "multiqc/fastq_trimmed/" + prefix +"_" + unique_id + ".html"
        list_inputs.append(directory_data_trimmed)
        list_inputs.append(html_trimmed)


use_mapping = config["USAGE"].get("MAPPING")
if use_mapping:
    STAR = config["DEPENDANCES"]["MAPPING"].get("STAR")
    if not STAR:
        STAR = "STAR"
    check_excutable(STAR)

    samtools = config["DEPENDANCES"]["MAPPING"].get("SAMTOOLS")
    if not samtools:
        samtools = "samtools"
    check_excutable("samtools")

    Rscript = config["DEPENDANCES"]["GENERAL"].get("RSCRIPT")
    if not Rscript:
        Rscript = "Rscript"
    check_excutable(Rscript)

    gff3 = config["GENERAL"].get("GFF3")
    if not os.path.isfile(gff3):
        gff3 = working_directory + "/1-raw_data/annotation/" + gff3

    if not os.path.isfile(gff3):
        print_once(f"Erreur : Le fichier gff3 : {gff3} n'existe pas.")
        sys.exit(1)
    else : 
        print_once(f"Fichier {gff3} ...... OK")

    mane = config["GENERAL"].get("MANE")
    if not os.path.isfile(mane):
        mane = working_directory + "/1-raw_data/annotation/" + mane

    if not os.path.isfile(mane):
        print_once(f"Erreur : Le fichier mane : {mane} n'existe pas.")
        sys.exit(1)
    else : 
        print_once(f"Fichier {mane} ...... OK")


    include: "modules/mapping.smk"
    print_once("Module mapping.smk ...... OK")

    if use_qc:
        bam_stat = config["DEPENDANCES"]["QC"].get("BAM_STAT")
        if not bam_stat:
            bam_stat = "bam_stat.py"
        check_excutable(bam_stat)

        read_GC = config["DEPENDANCES"]["QC"].get("READ_GC")
        if not read_GC:
            read_GC = "read_GC.py"
        check_excutable(read_GC)

        include: "modules/quality_control_bam.smk"
        print_once("Module quality_control_bam.smk ...... OK")

        directory_data_bam = path_qc + "multiqc/BAM/" + name_genome + "/" + prefix +"_" + unique_id + "_data/"
        html_bam = path_qc + "multiqc/BAM/" + name_genome + "/" + prefix +"_" + unique_id + ".html"
        list_inputs.append(directory_data_bam)
        list_inputs.append(html_bam)
    # list_inputs.append(count_report)
    # list_inputs.append(directory_count_results)

    #Vérification des dépendances SpliceLauncher
    use_SpliceLauncher = config["USAGE"].get("SPLICELAUNCHER")
    if use_SpliceLauncher:

        perl = config["DEPENDANCES"]["GENERAL"].get("PERL")
        if not perl:
            perl = "perl"
        check_excutable(perl)

        bedtools = config["DEPENDANCES"]["ANALYSES"].get("BEDTOOLS")
        if not bedtools:
            bedtools = "bedtools"
        check_excutable(bedtools)

        SpliceLauncher = config["DEPENDANCES"]["ANALYSES"].get("SPLICELAUNCHER")
        if not os.path.isdir(SpliceLauncher):
            print_once(f"Erreur : Le dossier {SpliceLauncher} n'existe pas.")
            sys.exit(1)
        else : 
            print_once(f"Dossier {SpliceLauncher} ...... OK")
        
        use_sashimi = config["SPLICELAUNCHER"]["SASHIMI_PLOT"].get("USE")
        if use_sashimi:
            python = config["DEPENDANCES"]["GENERAL"].get("PYTHON")
            if not python :
                python = "python"
            check_excutable(python)

            ggsashimi = config["DEPENDANCES"]["VISUALISATION"].get("GGSASHIMI")
            if not os.path.isfile(ggsashimi) :
                print_once(f"Erreur : path {ggsashimi} n'existe pas.")
                sys.exit(1)
            else : 
                print_once(f"Path {ggsashimi} existant ...... OK")


            directories_unique_junctions = expand(path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/sashimi_plot/unique_junctions/",reads= all_samples)
            directories_significant_junctions = expand(path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/sashimi_plot/significant_junctions/",reads= all_samples)
            list_inputs.append(directories_unique_junctions)
            list_inputs.append(directories_significant_junctions)

    
        include: "modules/SpliceLauncher.smk"
        print_once("Module SpliceLauncher ...... OK")

        count_report = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_report_" + date + ".txt"
        directory_count_results = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results"
        filterFilesSample = expand(path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.significant_junctions" + extension,reads= all_samples)
        filterFile = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + "_outputSpliceLauncher.significant_junctions" + extension,
        list_inputs.append(filterFile)
        list_inputs.append(filterFilesSample)

    use_variant_calling = config["USAGE"].get("VARIANT_CALLING")
    if use_variant_calling :
        gatk = config["DEPENDANCES"]["ANALYSES"].get("GATK")
        if not gatk:
            gatk = "gatk"
        check_excutable(gatk)

        bcftools = config["DEPENDANCES"]["GENERAL"].get("BCFTOOLS")
        if not bcftools:
            bcftools = "bcftools"
        check_excutable(bcftools)

        bgzip = config["DEPENDANCES"]["GENERAL"].get("BGZIP")
        if not bgzip:
            bgzip = "bgzip"
        check_excutable(bgzip)

        vcf_known_sites = config["GENERAL"].get("VCF_KNOWN_SITES")
        if not os.path.isfile(vcf_known_sites):
            vcf_known_sites = working_directory + "/1-raw_data/annotation/" + vcf_known_sites

        if not os.path.isfile(vcf_known_sites):
            print_once(f"Erreur : Le fichier vcf_known_sites : {vcf_known_sites} n'existe pas.")
            sys.exit(1)
        else : 
            print_once(f"Fichier {vcf_known_sites} ...... OK")

        include: "modules/gatk_recalibrated.smk"
        print_once("Module gatk_recalibrated ...... OK")
        
        include: "modules/gatk_variant_calling_germline.smk"
        print_once("Module gatk_variant_calling_germline ...... OK")

        # vcf = expand(path_results + "/GATK/variant_calling/{reads}.HaplotypeCaller.vcf.gz",reads= all_samples)
        # tbi = expand(path_results + "/GATK/variant_calling/{reads}.HaplotypeCaller.vcf.gz.tbi",reads= all_samples)
        # list_inputs.append(vcf)
        # list_inputs.append(tbi)
        vcf = expand(path_results + "/GATK/variant_calling/{reads}.filtered.vcf.gz",reads= all_samples)
        tbi = expand(path_results + "/GATK/variant_calling/{reads}.filtered.vcf.gz.tbi",reads= all_samples)
        list_inputs.append(vcf)
        list_inputs.append(tbi)

    if not list_inputs:
        bai = expand(path_bam + name_genome + "/mapping/{reads}.sorted.bam.bai",reads= all_samples)
        list_inputs.append(bai)

else:
    print_once("Nothing mapped ...... nothing analyzed")

if not list_inputs:
    read1 = expand(os.path.abspath(path_fastq + "{reads}.1.fastq.gz"),reads= all_samples)
    read2 = expand(os.path.abspath(path_fastq + "{reads}.2.fastq.gz"),reads= all_samples)
    list_inputs.append(read1)
    list_inputs.append(read2)

    if use_trimming:
        html = expand(path_qc + "fastp_trimming/{reads}.html",reads= all_samples)
        json = expand(path_qc + "fastp_trimming/{reads}.json",reads= all_samples)
        list_inputs.append(html)
        list_inputs.append(json)

os.environ["SNAKEMAKE_PRINT"] = "true"
#Calling Snakemake module
#include: "modules/featureCounts.smk"
#include: "modules/hisat2.smk"


# rule all:
#     input:
#         mergeFile = working_directory + prefix + "/2-Counts/merge.counts"

# rule all:
#     input:
#         # read1 =expand(os.path.abspath(working_directory + "/2-processed_data/samples/{reads}_1.fastq.gz"),reads= all_samples),
#         # read2 =expand(os.path.abspath(working_directory + "/2-processed_data/samples/{reads}_2.fastq.gz"),reads= all_samples),
#         #bam = expand(working_directory + prefix + "/1-mapping/{reads}.sorted.bam.bai",reads= all_samples),
#         directory_data = working_directory + "/2-processed_data/quality_control/multiqc/" + prefix +"_fastq_raw_" + unique_id + "_data/",
#         html =working_directory + "/2-processed_data/quality_control/multiqc/" + prefix  + "_fastq_raw_" + unique_id + ".html",
#         directory_data2 = working_directory + "/2-processed_data/quality_control/multiqc/" + prefix +"_fastq_trimmed_" + unique_id + "_data/",
#         html2 = working_directory  + "/2-processed_data/quality_control/multiqc/" + prefix  + "_fastq_trimmed_" + unique_id + ".html",
#         directory_data_bam = working_directory + prefix + "/quality_control/multiqc/" + prefix +"_bam_" + unique_id + "_data/",
#         html_bam = working_directory + prefix + "/quality_control/multiqc/" + prefix  + "_bam_" + unique_id + ".html",

rule all:
    input:
        list_inputs
        # bam = expand(working_directory + prefix + "/1-mapping/STAR/{reads}_ Aligned.out.bam",reads= all_samples),
        # outFileNamePrefix = expand(working_directory + prefix + "/1-mapping/STAR/{reads}_",reads= all_samples)
    shell:
        "bash scripts/clear_cache.sh"

# rule clear_cache:
#     output: touch(working_directory + "/log/clear_cache.done")
#     shell:
#         "bash scripts/clear_cache.sh"


#hook
onsuccess:
    memory_release()




onerror:
    mergeFile = path_results + "/SpliceLauncher/mergeFile/" + prefix + "_" + unique_id + ".txt "
    count_report = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "._report_" + date + ".txt "
    directory_count_results = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results "
    command =  "rm -rf " + mergeFile + count_report + directory_count_results + "&& bash scripts/clear_cache.sh"
    memory_release()
    shell(command)


# rule delete_tmp:
#     input: 
#         rules.all.output

#     params:
#         tmp_directory = working_directory + "/.tmp/"
    
#     run:
#         delete_directory(tmp_directory)