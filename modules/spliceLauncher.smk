rule spliceLauncherDB:
    input:
        gff3 =  gff3,
        mane = mane

    output:
        directory = directory(working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR"),
        bed = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/BEDannotation.bed",
        sjdb = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/SJDBannotation.sjdb",
        annot = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/SpliceLauncherAnnot.txt"

    params:
        spliceLauncher = spliceLauncher,
        Rscript = Rscript

    log:
        out = working_directory + "/log/spliceLauncherDB/"+ name_genome + ".out",
        err = working_directory + "/log/spliceLauncherDB/"+ name_genome + ".err"

    shell:
        "{params.Rscript} {params.spliceLauncher}/scripts/generateSpliceLauncherDB.r "
        "-i {input.gff3} "
        "-o {output.directory} "
        "--mane {input.mane} > {log.out} 2> {log.err}"



sjdbOverhang = config["SPLICELAUNCHER"]["INDEX"]["SJDB_OVERHANG"]
ram_byte = config["SPLICELAUNCHER"]["INDEX"]["RAM"]
genomeSAsparseD = config["SPLICELAUNCHER"]["INDEX"]["GENOME_SA_SPARSE_D"]
genomeSAindexNbases = config["SPLICELAUNCHER"]["INDEX"]["GENOME_SA_INDEX_NBASES"]
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
    reference = genome
    cmd = f"grep -v '>' {reference} | awk '{{ sum += length($0)}} END {{print sum}}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    genome_length = int(result.stdout.strip())

    # Calculer min(14, log2(genome_length) / 2 - 1)
    genomeSAindexNbases = min(14, math.log2(genome_length) / 2 - 1)

rule spliceLauncher_star_index:
    input:
        reference = genome,
        genomeDir = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR",
        sjdb = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/SJDBannotation.sjdb",
        gff3 = gff3
        

    output:
        reference = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/" + name_genome,
        SA = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/SA",
        SAindex = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/SAindex",
        chrLength = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/chrLength.txt",
        chrNameLength = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/chrNameLength.txt",
        chrName = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/chrName.txt",
        chrStart = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/chrStart.txt",
        exonGeTrInfo = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/exonGeTrInfo.tab",
        exonInfo = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/exonInfo.tab",
        geneInfo = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/geneInfo.tab",
        genome = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/Genome",
        genomeParameters = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/genomeParameters.txt",
        Log = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/Log.out",
        sjdbInfo = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/sjdbInfo.txt",
        sjdbListFromGTF = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/sjdbList.fromGTF.out.tab",
        sjdblist = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/sjdbList.out.tab",
        transcriptInfo = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/transcriptInfo.tab"

    params:
        star = STAR,
        sjdbOverhang = sjdbOverhang,
        genomeSAsparseD = genomeSAsparseD,
        limitGenomeGenerateRAM = ram_byte,
        genomeSAindexNbases = genomeSAindexNbases

    threads:
        config["SPLICELAUNCHER"]["INDEX"]["THREADS"]
    
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


def calculate_mem_per_rule(mem_mb, max_cores, threads):
    max_threads = max(int((int(max_cores)/ int(threads))),1)
    mem_per_process = int(mem_mb) / int(max_threads)
    return int(mem_per_process)


rule spliceLauncher_star_align:
    input:
        read1 = os.path.abspath(path_fastq + "{reads}.1.fastq.gz"),
        read2 = os.path.abspath(path_fastq + "{reads}.2.fastq.gz"),
        genomeDir = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR",
        SA = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/SA",
        SAindex = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/SAindex",

    output:
        bam = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}_Aligned.out.bam",
        log_final = path_bam + "spliceLauncher_STAR/" + name_genome + "/log_star/{reads}_Log.final.out",
        log = path_bam + "spliceLauncher_STAR/" + name_genome + "/log_star/{reads}_Log.out",
        log_progress = path_bam + "spliceLauncher_STAR/" + name_genome + "/log_star/{reads}_Log.progress.out",
        tab = path_bam + "spliceLauncher_STAR/" + name_genome + "/log_star/{reads}_SJ.out.tab",
        STARpass1 = directory(path_bam + "spliceLauncher_STAR/" + name_genome + "/log_star/{reads}__STARgenome"),
        STARgenome = directory(path_bam + "spliceLauncher_STAR/" + name_genome + "/log_star/{reads}__STARpass1")

    params:
        star = STAR,
        outFilterMismatchNmax = config["SPLICELAUNCHER"]["MAPPING"]["OUT_FILTER_MISMATCH_NMAX"],
        outFilterMultimapNmax = config["SPLICELAUNCHER"]["MAPPING"]["OUT_FILTER_MULTIMAP_NMAX"],
        outSJfilterIntronMaxVsReadN = config["SPLICELAUNCHER"]["MAPPING"]["OUT_SJ_FILTER_INTRON_MAXVSREADN"],
        outFileNamePrefix = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}_",
        outTmpDir = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}_tmp",
        directory_log = path_bam + "spliceLauncher_STAR/" + name_genome + "/log_star/",
        log_final = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}_Log.final.out",
        log = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}_Log.out",
        log_progress = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}_Log.progress.out",
        tab = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}_SJ.out.tab",
        STARpass1 = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}__STARgenome",
        STARgenome = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}__STARpass1"

    threads:
        config["SPLICELAUNCHER"]["MAPPING"]["THREADS"]
    
    resources:
        mem_mb = lambda wildcards,threads: calculate_mem_per_rule(mem_mb, max_cores, threads)

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
            "--genomeLoad NoSharedMemory")
            shell("mv {params.log} {params.log_final} {params.log_progress} {params.tab} {params.STARpass1} {params.STARgenome} {params.directory_log}")
        except Exception as e:
            shell("rm -rf {params.log} {params.log_final} {params.log_progress} {params.tab} {params.STARpass1} {params.STARgenome} {params.outTmpDir}")
            raise e


rule spliceLauncher_samtools_sort:
    input:
        bam = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}_Aligned.out.bam",

    output:
        sort_bam = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}.sorted.bam",
        csi = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}.sorted.bam.csi",

    params:
        samtools = samtools

    threads:
        config["PARAMS"]["STAR"]["THREADS"]

    shell:
        "{params.samtools} sort "
        "-o {output.sort_bam} "
        "--output-fmt BAM "
        "--write-index "
        "--threads {threads} "
        "{input.bam} && "
        "rm {input.bam}"

rule spliceLauncher_create_bed:
    input:
        bam = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}.sorted.bam",
        csi = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}.sorted.bam.csi",

    output:
        bed = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}_juncs.bed"

    params:
        samtools = samtools,
        bedtools = bedtools,
        spliceLauncher = spliceLauncher,
        perl = perl

    shell:
        "{params.samtools} view "
        "-b -f 0x40 -F 1024 {input.bam} | "
        "{params.bedtools} bamtobed -bed12 -i stdin | "
        "awk '{{if($10>1){{print $0}}}}' | "
        "{params.perl} {params.spliceLauncher}/scripts/bedBlocks2IntronsCoords.pl y - | "
        "awk '{{if($5==255){{print $0}}}}' > {output.bed} && "
        "{params.samtools} view "
        "-b -f 0x80 -F 1024 {input.bam} | "
        "{params.bedtools} bamtobed -bed12 -i stdin | "
        "awk '{{if($10>1){{print $0}}}}' | "
        "{params.perl} {params.spliceLauncher}/scripts/bedBlocks2IntronsCoords.pl n - | "
        "awk '{{if($5==255){{print $0}}}}' >> {output.bed}"


rule spliceLauncher_count_Junctions:
    input:
        bed_ref = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/BEDannotation.bed",
        bed = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}_juncs.bed"

    output:
        counts = path_bam + "spliceLauncher_STAR/" + name_genome + "/getClosestExons/{reads}.count"

    params:
        bedtools = bedtools

    shell:
        "sort -k1,1 -k2,2n {input.bed} | "
        "uniq -c | "
        "awk 'BEGIN{{OFS=\"\\t\"}}{{print $2,$3,$4,$1,$6,$7}}' | "
        "{params.bedtools} intersect -s -wa -wb -a stdin -b {input.bed_ref} | "
        "awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,$3,$6,$10,$4}}' "
        "> {output.counts}"


genes_of_interest = config["SPLICELAUNCHER"]["ANALYSE"]["GENES_OF_INTEREST"]
if not os.path.isfile(genes_of_interest):
    genes_of_interest = working_directory + "/metadata/genes_of_interest/" + genes_of_interest

if not os.path.isfile(genes_of_interest):
    genes_of_interest = ""

else :
    print_once(f"Fichier {genes_of_interest} ...... OK")
    genes_of_interest = "-g " + genes_of_interest


rule spliceLauncher_merge_count:
    input:
        counts = expand(path_bam + "spliceLauncher_STAR/" + name_genome + "/getClosestExons/{reads}.count",reads= all_samples)

    output:
        mergeFile = path_results + "/spliceLauncher/mergeFile/" + prefix + "_" + unique_id + ".txt"
    
    params:
        spliceLauncher = spliceLauncher,
        perl = perl,
        directory = directory(path_results + "/spliceLauncher/mergeFiles/"),
        genes_of_interest = genes_of_interest

    shell:
        "{params.perl} {params.spliceLauncher}/scripts/joinJuncFiles.pl -c {input.counts} {params.genes_of_interest} > {output.mergeFile}"


sampleNames = result = '|'.join(all_samples)

nbIntervals = config["SPLICELAUNCHER"]["ANALYSE"]["NB_INTERVALS"]
if not nbIntervals:
    nbIntervals = int(10)

threshold = config["SPLICELAUNCHER"]["ANALYSE"]["THRESHOLD"]
if not threshold:
    threshold = int(1)

min_cov = config["SPLICELAUNCHER"]["ANALYSE"]["MIN_COV"]
if not min_cov:
    min_cov = int(5)

transcriptList = config["SPLICELAUNCHER"]["ANALYSE"]["TRANSCRIPT_LIST"]
removeOther = config["SPLICELAUNCHER"]["ANALYSE"]["REMOVE_OTHER"]
if transcriptList and os.path.isfile(transcriptList):
    transcriptList = "--TranscriptList " + transcriptList + " "
    if removeOther:
        removeOther = "--removeOther "
    else:
        removeOther = ""
else:
    transcriptList = ""
    removeOther = ""

txt = config["SPLICELAUNCHER"]["ANALYSE"]["TEXT"]
if txt:
    txt = "--text "
    extension = ".txt"

else:
    txt = ""
    extension = ".xslx"

bedOut = config["SPLICELAUNCHER"]["ANALYSE"]["BED_OUT"]
if bedOut:
    bedOut = "--bedOut "
else:
    bedOut = ""

Graphics =config["SPLICELAUNCHER"]["ANALYSE"]["GRAPHICS"]
if Graphics:
    Graphics = "--Graphics "

current_time = time.localtime()
#Formater la date et l'heure comme identifiant unique
date = time.strftime("%m-%d-%Y", current_time)


rule spliceLauncher_Analyse:
    input:
        mergeFile = path_results + "/spliceLauncher/mergeFile/" + prefix + "_" + unique_id + ".txt",
        annot = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/SpliceLauncherAnnot.txt"

    output:
        count_report = path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_report_" + date + ".txt",
        outputSpliceLauncher = path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + "_outputSpliceLauncher" + extension,
        directy_count_results = directory(path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results")

    params:
        spliceLauncher = spliceLauncher,
        Rscript = Rscript,
        nbIntervals = nbIntervals,
        threshold = threshold,
        min_cov = min_cov,
        transcriptList = transcriptList,
        removeOther = removeOther,
        txt = txt,
        bedOut = bedOut,
        Graphics = Graphics,
        directory = path_results + "/spliceLauncher/",
        sampleNames = sampleNames


    shell:
        "{params.Rscript} {params.spliceLauncher}/scripts/SpliceLauncherAnalyse.r "
        "--input {input.mergeFile} "
        "-O {params.directory} "
        "--RefSeqAnnot {input.annot} "
        "-n {params.nbIntervals} "
        "{params.transcriptList}"
        "{params.removeOther}"
        "{params.txt}"
        "{params.bedOut}"
        "{params.Graphics}"
        "--SampleNames '{params.sampleNames}' "
        "--threshold {params.threshold} "
        "--min_cov {params.min_cov} "


# length_all_samples = len(all_samples)

# rule spliceLauncher_filter_analyse:
#     input:
#         directy_count_results = path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results",
#         outputSpliceLauncher = path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + "_outputSpliceLauncher" + extension,

#     output:
#         filterFile = path_results + "/spliceLauncher/mergeFile/" + prefix + "_" + unique_id + ".filter" + extension
    
#     params:
#         Rscript = Rscript,
#         length_all_samples = length_all_samples,
#         txt = txt,
#         Graphics = Graphics

#     shell:
#         "{params.Rscript} scripts/spliceLaucher_filter_analyse.r "
#         "--input {input.outputSpliceLauncher} "
#         "--length {params.length_all_samples} "
#         "{params.txt}"
#         "{params.Graphics}"        

#create_directory_if_not_exists(params["directory"])
