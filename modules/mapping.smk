rule SpliceLauncherDB:
    input:
        gff3 =  gff3,
        mane = mane

    output:
        directory = directory(working_directory + "/2-processed_data/reference/" + name_genome),
        bed = working_directory + "/2-processed_data/reference/" + name_genome + "/BEDannotation.bed",
        sjdb = working_directory + "/2-processed_data/reference/" + name_genome + "/SJDBannotation.sjdb",
        annot = working_directory + "/2-processed_data/reference/" + name_genome + "/SpliceLauncherAnnot.txt"

    params:
        Rscript = Rscript

    log:
        out = working_directory + "/log/SpliceLauncherDB/"+ name_genome + ".out",
        err = working_directory + "/log/SpliceLauncherDB/"+ name_genome + ".err"

    shell:
        "{params.Rscript} scripts/generateSpliceLauncherDB.r "
        "-i {input.gff3} "
        "-o {output.directory} "
        "--mane {input.mane} > {log.out} 2> {log.err}"



def estimate_memory(genome_length, bytes_per_base = 10):
    memory_required_bytes = genome_length * bytes_per_base

    memory_required_mb = memory_required_bytes / (1024 * 1024)
    return memory_required_mb

sjdbOverhang = config["MAPPING"]["INDEX"].get("SJDB_OVERHANG")
ram_byte = config["MAPPING"]["INDEX"].get("RAM")
genomeSAsparseD = config["MAPPING"]["INDEX"].get("GENOME_SA_SPARSE_D")
genomeSAindexNbases = config["MAPPING"]["INDEX"].get("GENOME_SA_INDEX_NBASES")
if not sjdbOverhang:
    sjdbOverhang = int(99)

if not ram_byte:
    ram_byte = int(total_ram * 0.5)
    mem_mb = int(ram_byte / (1024 * 1024))
else :
    mem_mb = int(int(ram_byte) / (1024 * 1024))

if not genomeSAsparseD:
    genomeSAsparseD=int(31000000000 / ram_byte) + 1


genome_length = calculate_genome_length(genome)
memory_needed_mb = estimate_memory(genome_length, 10)
if not genomeSAindexNbases:
    # Calculer min(14, log2(genome_length) / 2 - 1)
    genomeSAindexNbases = min(14, math.log2(genome_length) / 2 - 1)


rule star_index:
    input:
        reference = genome,
        genomeDir = working_directory + "/2-processed_data/reference/" + name_genome,
        sjdb = working_directory + "/2-processed_data/reference/" + name_genome + "/SJDBannotation.sjdb",
        gff3 = gff3
        

    output:
        reference = working_directory + "/2-processed_data/reference/" + name_genome + "/" + name_genome,
        SA = working_directory + "/2-processed_data/reference/" + name_genome + "/SA",
        SAindex = working_directory + "/2-processed_data/reference/" + name_genome + "/SAindex",
        chrLength = working_directory + "/2-processed_data/reference/" + name_genome + "/chrLength.txt",
        chrNameLength = working_directory + "/2-processed_data/reference/" + name_genome + "/chrNameLength.txt",
        chrName = working_directory + "/2-processed_data/reference/" + name_genome + "/chrName.txt",
        chrStart = working_directory + "/2-processed_data/reference/" + name_genome + "/chrStart.txt",
        exonGeTrInfo = working_directory + "/2-processed_data/reference/" + name_genome + "/exonGeTrInfo.tab",
        exonInfo = working_directory + "/2-processed_data/reference/" + name_genome + "/exonInfo.tab",
        geneInfo = working_directory + "/2-processed_data/reference/" + name_genome + "/geneInfo.tab",
        genome = working_directory + "/2-processed_data/reference/" + name_genome + "/Genome",
        genomeParameters = working_directory + "/2-processed_data/reference/" + name_genome + "/genomeParameters.txt",
        Log = working_directory + "/2-processed_data/reference/" + name_genome + "/Log.out",
        sjdbInfo = working_directory + "/2-processed_data/reference/" + name_genome + "/sjdbInfo.txt",
        sjdbListFromGTF = working_directory + "/2-processed_data/reference/" + name_genome + "/sjdbList.fromGTF.out.tab",
        sjdblist = working_directory + "/2-processed_data/reference/" + name_genome + "/sjdbList.out.tab",
        transcriptInfo = working_directory + "/2-processed_data/reference/" + name_genome + "/transcriptInfo.tab"

    params:
        star = STAR,
        sjdbOverhang = sjdbOverhang,
        genomeSAsparseD = genomeSAsparseD,
        limitGenomeGenerateRAM = ram_byte,
        genomeSAindexNbases = genomeSAindexNbases

    threads:
        config["MAPPING"]["INDEX"]["THREADS"]
    
    resources:
        mem_mb = mem_mb

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


def calculate_optimal_threads(mem_mb, memory_needed_mb, max_cores, requested_threads):
    total_mem_needed = memory_needed_mb * (max_cores // requested_threads)
    if total_mem_needed > mem_mb :
        max_parallel_rules = mem_mb // memory_needed_mb
        optimal_threads = max(1, int(max_cores / max_parallel_rules))
    else:
        optimal_threads = requested_threads
    
    optimal_threads = min(optimal_threads, max_cores)
    return optimal_threads

def calculate_mem_per_rule(mem_mb, memory_needed_mb, optimal_threads):
    max_parallel_rules = mem_mb // memory_needed_mb
    mem_per_process = int(mem_mb / max_parallel_rules)
    return int(mem_per_process)

outFilterMismatchNmax = config["MAPPING"]["ALIGN"].get("OUT_FILTER_MISMATCH_NMAX")
if not outFilterMismatchNmax:
    outFilterMismatchNmax = int(2)

outFilterMultimapNmax = config["MAPPING"]["ALIGN"].get("OUT_FILTER_MULTIMAP_NMAX")
if not outFilterMultimapNmax:
    outFilterMultimapNmax = int(1)

outSJfilterIntronMaxVsReadN = config["MAPPING"]["ALIGN"].get("OUT_SJ_FILTER_INTRON_MAXVSREADN")
if not outSJfilterIntronMaxVsReadN:
    outSJfilterIntronMaxVsReadN = int(500000)

rule star_align:
    input:
        read1 = os.path.abspath(path_fastq + "{reads}.1.fastq.gz"),
        read2 = os.path.abspath(path_fastq + "{reads}.2.fastq.gz"),
        genomeDir = working_directory + "/2-processed_data/reference/" + name_genome,
        SA = working_directory + "/2-processed_data/reference/" + name_genome + "/SA",
        SAindex = working_directory + "/2-processed_data/reference/" + name_genome + "/SAindex",

    output:
        sort_bam = os.path.abspath(path_bam + name_genome + "/mapping/{reads}.sorted.bam"),
        bai = path_bam + name_genome + "/mapping/{reads}.sorted.bam.bai",
        csi = path_bam + name_genome + "/mapping/{reads}.sorted.bam.csi",
        log_final = path_bam  + name_genome + "/mapping/log_star/{reads}_Log.final.out",
        log = path_bam  + name_genome + "/mapping/log_star/{reads}_Log.out",
        log_progress = path_bam  + name_genome + "/mapping/log_star/{reads}_Log.progress.out",
        tab = path_bam  + name_genome + "/mapping/log_star/{reads}_SJ.out.tab",
        STARpass1 = directory(path_bam  + name_genome + "/mapping/log_star/{reads}__STARgenome"),
        STARgenome = directory(path_bam  + name_genome + "/mapping/log_star/{reads}__STARpass1")

    params:
        star = STAR,
        samtools = samtools,
        bam = os.path.abspath(path_bam  + name_genome + "/mapping/{reads}_Aligned.out.bam"),
        outFilterMismatchNmax = outFilterMismatchNmax,
        outFilterMultimapNmax = outFilterMultimapNmax,
        outSJfilterIntronMaxVsReadN = outFilterMultimapNmax,
        outFileNamePrefix = path_bam  + name_genome + "/mapping/{reads}_",
        outTmpDir = path_bam  + name_genome + "/mapping/{reads}_tmp",
        directory_log = path_bam  + name_genome + "/mapping/log_star/",
        log_final = path_bam  + name_genome + "/mapping/{reads}_Log.final.out",
        log = path_bam  + name_genome + "/mapping/{reads}_Log.out",
        log_progress = path_bam  + name_genome + "/mapping/{reads}_Log.progress.out",
        tab = path_bam  + name_genome + "/mapping/{reads}_SJ.out.tab",
        STARpass1 = path_bam  + name_genome + "/mapping/{reads}__STARgenome",
        STARgenome = path_bam  + name_genome + "/mapping/{reads}__STARpass1",

    threads:
        lambda wildcards: calculate_optimal_threads(mem_total_mb, memory_needed_mb, max_cores, int(config["MAPPING"]["ALIGN"]["THREADS"]))
    
    resources:
        mem_mb = lambda wildcards,threads: calculate_mem_per_rule(mem_total_mb, memory_needed_mb, threads)

    run:
        create_directory_if_not_exists(params["directory_log"])
        try:
            shell("{params.star} "
            "--runThreadN {threads} "
            "--outSAMstrandField intronMotif "
            "--outFilterMismatchNmax {params.outFilterMismatchNmax} "
            "--outFilterMultimapNmax {params.outFilterMultimapNmax} "
            "--genomeDir {input.genomeDir} "
            "--readFilesIn {input.read1} {input.read2} "
            "--readFilesCommand zcat "
            "--outSAMunmapped Within "
            "--outSAMtype BAM Unsorted "
            "--outSJfilterOverhangMin -1 8 8 8 "
            "--outSJfilterCountUniqueMin -1 1 1 1 "
            "--outSJfilterDistToOtherSJmin 0 0 0 0 "
            "--alignSJstitchMismatchNmax 0 -1 -1 -1 "
            "--outSJfilterIntronMaxVsReadN {params.outSJfilterIntronMaxVsReadN} "
            "--twopassMode Basic "
            "--outTmpDir {params.outTmpDir} "
            "--outSAMheaderHD \\@HD VN:1.4 SO:Unsorted "
            "--outFileNamePrefix {params.outFileNamePrefix} "
            "--genomeLoad NoSharedMemory ")
            shell("{params.samtools} sort -@ {threads} -o {output.sort_bam} --output-fmt BAM --write-index {params.bam}")
            shell("{params.samtools} index -@ {threads} -b {output.sort_bam}")
            shell("mv {params.log} {params.log_final} {params.log_progress} {params.tab} {params.STARpass1} {params.STARgenome} {params.directory_log}")
            shell("rm -rf {params.bam}")
        except Exception as e:
            shell("rm -rf {params.log} {params.log_final} {params.log_progress} {params.tab} {params.STARpass1} {params.STARgenome} {params.outTmpDir} {params.bam}")
            raise e

# rule samtools_sort_bam:
#     input:
#         bam = os.path.abspath(path_bam  + name_genome + "/{reads}_Aligned.out.bam"),

#     output:
#         sort_bam = os.path.abspath(path_bam + name_genome + "/{reads}.sorted.bam"),
#         csi = path_bam + name_genome + "/{reads}.sorted.bam.csi",

#     params:
#         samtools = samtools

#     threads:
#         config["MAPPING"]["ALIGN"]["THREADS"]

#     shell:
#         "{params.samtools} sort "
#         "-o {output.sort_bam} "
#         "--output-fmt BAM "
#         "--write-index "
#         "--threads {threads} "
#         "{input.bam} "

# rule samtools_index_bam:
#     input:
#         sort_bam = os.path.abspath(path_bam + name_genome + "/{reads}.sorted.bam")

#     output:
#         bai = path_bam + name_genome + "/{reads}.sorted.bam.bai",

#     params:
#         samtools = samtools

#     threads:
#         config["MAPPING"]["ALIGN"]["THREADS"]

#     shell:
#         "{params.samtools} index "
#         "-@ {threads}"
#         "-b {input.sort_bam}"