rule spliceLaucherDB:
    input:
        gff3 =  config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/annotation/" + config["GENERAL"]["DATA_INPUTS"]["GFF3"],
        mane = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/annotation/" + config["GENERAL"]["DATA_INPUTS"]["MANE"]

    output:
        directory = directory(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR"),
        bed = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/BEDannotation.bed",
        sjdb = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/SJDBannotation.sjdb",
        annot = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/SpliceLauncherAnnot.txt"

    params:
        spliceLaucher = config["DEPENDANCES"]["ANALYSES"]["SPLICELAUCHER"]

    log:
        out = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/spliceLaucherDB/"+ config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] +".out",
        err = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/spliceLaucherDB/"+ config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] +".err"

    shell:
        "Rscript {params.spliceLaucher}/scripts/generateSpliceLauncherDB.r "
        "-i {input.gff3} "
        "-o {output.directory} "
        "--mane {input.mane} > {log.out} 2> {log.err}"



sjdbOverhang = config["SPLICELAUCHER"]["INDEX"]["SJDBOVERHANG"]
ram_byte = config["SPLICELAUCHER"]["INDEX"]["RAM"]
genomeSAsparseD = config["SPLICELAUCHER"]["INDEX"]["GENOMESASPARSED"]
genomeSAindexNbases = config["SPLICELAUCHER"]["INDEX"]["GENOMESAINDEXNBASES"]
if not sjdbOverhang:
    sjdbOverhang = int(99)

if not ram_byte:
    cmd = "free -tb | grep 'Mem:' | awk ' { print $2 } ' | awk '{printf \"%.f\", $1*0.5}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    ram_byte = int(result.stdout.strip())
    cmd = "free -tm | grep 'Mem:' | awk ' { print $2 } ' | awk '{printf \"%.f\", $1*0.5}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    mem_mb = int(result.stdout.strip())
else :
    mem_mb = int((int(ram_byte) / 1024))

if not genomeSAsparseD:
    genomeSAsparseD=int(31000000000 / ram_byte) + 1

if not genomeSAindexNbases:
    reference = os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"])
    cmd = f"grep -v '>' {reference} | awk '{{ sum += length($0)}} END {{print sum}}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    genome_length = int(result.stdout.strip())

    # Calculer min(14, log2(genome_length) / 2 - 1)
    genomeSAindexNbases = min(14, math.log2(genome_length) / 2 - 1)

rule star_index:
    input:
        reference = os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"]),
        genomeDir = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR",
        sjdb = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/SJDBannotation.sjdb",
        gff3 =  config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/annotation/" + config["GENERAL"]["DATA_INPUTS"]["GFF3"]
        

    output:
        reference = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"],
        SA = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/SA",
        SAindex = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/SAindex",
        chrLength = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/chrLength.txt",
        chrNameLength = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/chrNameLength.txt",
        chrName = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/chrName.txt",
        chrStart = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/chrStart.txt",
        exonGeTrInfo = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/exonGeTrInfo.tab",
        exonInfo = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/exonInfo.tab",
        geneInfo = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/geneInfo.tab",
        genome = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/Genome",
        genomeParameters = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/genomeParameters.txt",
        Log = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/Log.out",
        sjdbInfo = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/sjdbInfo.txt",
        sjdbListFromGTF = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/sjdbList.fromGTF.out.tab",
        sjdblist = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/sjdbList.out.tab",
        transcriptInfo = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/transcriptInfo.tab"

    params:
        star = config["DEPENDANCES"]["MAPPING"]["STAR"],
        reference = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "_STAR/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"],
        sjdbOverhang = sjdbOverhang,
        genomeSAsparseD = genomeSAsparseD,
        limitGenomeGenerateRAM = ram_byte,
        genomeSAindexNbases = genomeSAindexNbases

    threads:
        config["PARAMS"]["STAR"]["THREADS"]
    
    resources:
        mem_mb = mem_mb

    log:

    shell:
        "ln -sfn {input.reference} {output.reference} && "
        "{params.star} "
        "--runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir {input.genomeDir} "
        "--genomeFastaFiles {output.reference} "
        "--sjdbFileChrStartEnd {input.sjdb} "
        "--sjdbGTFfile {input.gff3} "
        "--sjdbOverhang {params.sjdbOverhang} "
        "--genomeSAsparseD {params.genomeSAsparseD} "
        "--limitGenomeGenerateRAM {params.limitGenomeGenerateRAM} "
        "--genomeSAindexNbases {params.genomeSAindexNbases}"



# rule star_align:
#     input:
#         trimmed_read1 = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_trimmed/{reads}_1.trimmed.fastq.gz",
#         trimmed_read2 = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_trimmed/{reads}_2.trimmed.fastq.gz",
#         genomeDir = 

#     output:

#     params:
#         star = config["DEPENDANCES"]["MAPPING"]["STAR"],
#         outFilterMismatchNmax =
#         outFilterMultimapNmax = 
#         outTmpDir = 
#         outFileNamePrefix = 
#         outSJfilterIntronMaxVsReadN =


#     threads:
#         config["PARAMS"]["STAR"]["THREADS"]


#     log:

#     run:
#     "{params.star} "
#     "--runThreadN {threads} "
#     "--outSAMstrandField intronMotif "
#     "--outFilterMismatchNmax {params.outFilterMismatchNmax} "
#     "--outFilterMultimapNmax {params.outFilterMultimapNmax} "
#     "--genomeDir {input.genomeDir} "
#     "--readFilesIn {input.trimmed_read1} {input.trimmed_read2} "
#     "--readFilesCommand zcat "
#     "--outSAMunmapped Within "
#     "--outSAMtype BAM Unsorted "
#     "--outSJfilterOverhangMin -1 8 8 8 "
#     "--outSJfilterCountUniqueMin -1 1 1 1 "
#     "--outSJfilterDistToOtherSJmin 0 0 0 0 "
#     "--alignSJstitchMismatchNmax 0 -1 -1 -1 "
#     "--outSJfilterIntronMaxVsReadN {params.outSJfilterIntronMaxVsReadN} "
#     "--twopassMode Basic "
#     "--outTmpDir {params.outTmpDir} "
#     "--outSAMheaderHD \@HD VN:1.4 SO:Unsorted "
#     "--outFileNamePrefix {params.outFileNamePrefix} "
#     "--genomeLoad NoSharedMemory "
    # """STAR \
    # --outSAMstrandField intronMotif \
    # --outFilterMismatchNmax 2 \
    # --outFilterMultimapNmax 10 \
    # --genomeDir genome/GRCh37 \
    # --readFilesIn /home/micro/Marco_Corentin/06-data/02-fastq/01-ROSALIND/run_fastq_2023_04/1409241167_S6_L001_R1_001B.fastq.gz /home/micro/Marco_Corentin/06-data/02-fastq/01-ROSALIND/run_fastq_2023_04/1409241167_S6_L001_R2_001A.fastq.gz \
    # --readFilesCommand zcat \
    # --runThreadN 3 \
    # --outSAMunmapped Within \
    # --outSAMtype BAM Unsorted \
    # --outSJfilterOverhangMin -1 8 8 8 \
    # --outSJfilterCountUniqueMin -1 1 1 1 \
    # --outSJfilterCountTotalMin -1 1 1 1 \
    # --outSJfilterDistToOtherSJmin 0 0 0 0 \
    # --alignSJstitchMismatchNmax 0 -1 -1 -1 \
    # --outSJfilterIntronMaxVsReadN 500000 \
    # --twopassMode Basic \
    # --outTmpDir toto \
    # --outSAMheaderHD \@HD VN:1.4 SO:Unsorted \
    # --outFileNamePrefix test_bam \
    # --genomeLoad NoSharedMemory > test_bam.log 2>&1"""



# rule samtools_index:


# rule create_bed:


# rule count_Junctions:


# rule merge_count:

