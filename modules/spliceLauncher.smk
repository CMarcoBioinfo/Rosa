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



def estimate_star_memory(genome_length):
    bytes_per_base = 10
    memory_required_bytes = genome_length * bytes_per_base

    memory_required_mb = memory_required_bytes / (1024 * 1024)
    return memory_required_mb

sjdbOverhang = config["SPLICELAUNCHER"]["INDEX"].get("SJDB_OVERHANG")
ram_byte = config["SPLICELAUNCHER"]["INDEX"].get("RAM")
genomeSAsparseD = config["SPLICELAUNCHER"]["INDEX"].get("GENOME_SA_SPARSE_D")
genomeSAindexNbases = config["SPLICELAUNCHER"]["INDEX"].get("GENOME_SA_INDEX_NBASES")
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
memory_needed_mb = estimate_star_memory(genome_length)
if not genomeSAindexNbases:
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

def calculate_mem_per_rule(mem_mb, optimal_threads):
    max_parallel_rules = mem_mb // memory_needed_mb
    mem_per_process = int(mem_mb / max_parallel_rules)
    return int(mem_per_process)

outFilterMismatchNmax = config["SPLICELAUNCHER"]["MAPPING"].get("OUT_FILTER_MISMATCH_NMAX")
if not outFilterMismatchNmax:
    outFilterMismatchNmax = int(2)

outFilterMultimapNmax = config["SPLICELAUNCHER"]["MAPPING"].get("OUT_FILTER_MULTIMAP_NMAX")
if not outFilterMultimapNmax:
    outFilterMultimapNmax = int(1)

outSJfilterIntronMaxVsReadN = config["SPLICELAUNCHER"]["MAPPING"].get("OUT_SJ_FILTER_INTRON_MAXVSREADN")
if not outSJfilterIntronMaxVsReadN:
    outSJfilterIntronMaxVsReadN = int(500000)


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
        outFilterMismatchNmax = outFilterMismatchNmax,
        outFilterMultimapNmax = outFilterMultimapNmax,
        outSJfilterIntronMaxVsReadN = outFilterMultimapNmax,
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
        lambda wildcards: calculate_optimal_threads(mem_total_mb, memory_needed_mb, max_cores, int(config["SPLICELAUNCHER"]["MAPPING"]["THREADS"]))
    
    resources:
        mem_mb = lambda wildcards,threads: calculate_mem_per_rule(mem_total_mb,threads)

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
        sort_bam = os.path.abspath(path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}.sorted.bam"),
        csi = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}.sorted.bam.csi",

    params:
        samtools = samtools

    threads:
        config["SPLICELAUNCHER"]["MAPPING"]["THREADS"]

    shell:
        "{params.samtools} sort "
        "-o {output.sort_bam} "
        "--output-fmt BAM "
        "--write-index "
        "--threads {threads} "
        "{input.bam} && "
        "rm {input.bam}"

rule spliceLauncher_samtools_index:
    input:
        bam = os.path.abspath(path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}.sorted.bam"),

    output:
        bai = path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}.sorted.bam.bai",

    params:
        samtools = samtools

    shell:
        "{params.samtools} index "
        "-b {input.bam}"


rule spliceLauncher_create_bed:
    input:
        bam = os.path.abspath(path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}.sorted.bam"),
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


genes_of_interest = config["SPLICELAUNCHER"]["ANALYSE"].get("GENES_OF_INTEREST")
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
nbIntervals = config["SPLICELAUNCHER"]["ANALYSE"].get("NB_INTERVALS")
if not nbIntervals:
    nbIntervals = int(10)

threshold = config["SPLICELAUNCHER"]["ANALYSE"].get("THRESHOLD")
if not threshold:
    threshold = int(1)

min_cov = config["SPLICELAUNCHER"]["ANALYSE"].get("MIN_COV")
if not min_cov:
    min_cov = int(5)

transcriptList = config["SPLICELAUNCHER"]["ANALYSE"].get("TRANSCRIPT_LIST")
removeOther = config["SPLICELAUNCHER"]["ANALYSE"].get("REMOVE_OTHER")
if transcriptList and os.path.isfile(transcriptList):
    transcriptList = "--TranscriptList " + transcriptList + " "
    if removeOther:
        removeOther = "--removeOther "
    else:
        removeOther = ""
else:
    transcriptList = ""
    removeOther = ""

txt = config["SPLICELAUNCHER"]["ANALYSE"].get("TXT")
if txt:
    txt = "--text "
    extension = ".txt"

else:
    txt = ""
    extension = ".xlsx"

#Formater la date et l'heure comme identifiant unique
current_time = time.localtime()
date = time.strftime("%m-%d-%Y", current_time)
listOutputSpliceLauncher = []
count_report = path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_report_" + date + ".txt",
listOutputSpliceLauncher.append(count_report)
outputSpliceLauncher = path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + "_outputSpliceLauncher" + extension,
listOutputSpliceLauncher.append(outputSpliceLauncher)

bedOut = config["SPLICELAUNCHER"]["ANALYSE"].get("BED_OUT")
if bedOut:
    bedOut = "--bedOut "
    bed = path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + ".bed"
    listOutputSpliceLauncher.append(bed)

else:
    bedOut = ""

Graphics = config["SPLICELAUNCHER"]["ANALYSE"].get("GRAPHICS")
if Graphics:
    Graphics = "--Graphics "
    pdf = expand(path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.pdf",reads= all_samples)
    pdf_significant_genes = expand(path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.significant_genes.pdf",reads= all_samples)
    pdf_significant_junctions = expand(path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.significant_junctions.pdf",reads= all_samples)
    listOutputSpliceLauncher.append(pdf)
    listOutputSpliceLauncher.append(pdf_significant_genes)
    listOutputSpliceLauncher.append(pdf_significant_junctions)




rule spliceLauncher_Analyse:
    input:
        mergeFile = path_results + "/spliceLauncher/mergeFile/" + prefix + "_" + unique_id + ".txt",
        annot = working_directory + "/2-processed_data/reference/" + name_genome + "_spliceLauncher_STAR/SpliceLauncherAnnot.txt"

    output:
        listOutputSpliceLauncher

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
        "--min_cov {params.min_cov}"


length_all_samples = len(all_samples)
outputFilterAnalyse = []
filterFilesSampleSignificant = expand(path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.significant_junctions" + extension,reads= all_samples)

filterFileSignificant = path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + "_outputSpliceLauncher.significant_junctions" + extension,
filterFileUnique = path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + "_outputSpliceLauncher.unique_junctions" + extension,
outputFilterAnalyse.append(filterFileSignificant)
outputFilterAnalyse.append(filterFileUnique)
outputFilterAnalyse.append(filterFilesSampleSignificant)

filterFilesSampleSignificantFilter = expand(path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.significant_junctions.filter" + extension,reads= all_samples)
filterFilesSampleUniqueFilter = expand(path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.unique_junctions.filter" + extension,reads= all_samples)
outputFilterAnalyse.append(filterFilesSampleSignificantFilter)
outputFilterAnalyse.append(filterFilesSampleUniqueFilter)

maxSignificantSamples = config["SPLICELAUNCHER"]["POST_ANALYSE"]["SIGNIGICANT_JUNCTIONS"].get("MAX_SAMPLES")
if not maxSignificantSamples :
    maxSignificantSamples = int(2)

# minSignificantReads = config["SPLICELAUNCHER"]["POST_ANALYSE"]["SIGNIGICANT_JUNCTIONS"]["MIN_READS_SAMPLE"]
# if not minSignificantReads :
#     minSignificantReads = int(10)

maxUniqueSamples = config["SPLICELAUNCHER"]["POST_ANALYSE"]["UNIQUE_JUNCTIONS"].get("MAX_SAMPLES")
if not maxUniqueSamples : 
    maxUniqueSamples = int(2)

minUniqueReads = config["SPLICELAUNCHER"]["POST_ANALYSE"]["UNIQUE_JUNCTIONS"].get("MIN_READS_SAMPLE")
if not minUniqueReads :
    minUniqueReads = int(10)


rule spliceLauncher_filter_analyse:
    input:
        outputSpliceLauncher = path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + "_outputSpliceLauncher" + extension,

    output:
        outputFilterAnalyse

    params:
        Rscript = Rscript,
        length_all_samples = length_all_samples,
        txt = txt,
        Graphics = Graphics,
        sampleNames = sampleNames,
        maxSignificantSamples = maxSignificantSamples,
        #minSignificantReads = minSignificantReads,
        maxUniqueSamples = maxUniqueSamples,
        minUniqueReads = minUniqueReads,
        directory = path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/"

    shell:
        "{params.Rscript} scripts/spliceLauncher_filter_analyse.r "
        "--input {input.outputSpliceLauncher} "
        "--outputSignificant {output[0]} "
        "--outputUnique {output[1]} "
        "--length {params.length_all_samples} "
        #"--minSignificantReads {params.minSignificantReads} "
        "--maxSignificantSamples {params.maxSignificantSamples} "
        "--minUniqueReads {params.minUniqueReads} "
        "--maxUniqueSamples {params.maxUniqueSamples} "
        "--directory {params.directory} "
        "{params.txt}"
        "{params.Graphics}"
        "--SampleNames '{params.sampleNames}' "  


gtf = config["GENERAL"].get("GTF")
if not os.path.isfile(gtf):
    gtf = working_directory + "/1-raw_data/annotation/" + gtf
if not os.path.isfile(gtf):
    gtf = -1

color = config["SPLICELAUNCHER"]["SASHIMI_PLOT"].get("COLOR")
if not os.path.isfile(color):
    gtf = working_directory + "/metadata/other/" + color
if not os.path.isfile(color):
    color = -1

extend_bp = config["SPLICELAUNCHER"]["SASHIMI_PLOT"].get("EXTEND_BP")
if not extend_bp:
    extend_bp = 50

MinThresholdNbReads = config["SPLICELAUNCHER"]["SASHIMI_PLOT"].get("MIN_THRESHOLD_NB_READS")
if not MinThresholdNbReads:
    MinThresholdNbReads = -1

nb_samples = config["SPLICELAUNCHER"]["SASHIMI_PLOT"].get("NUMBER_SAMPLES")
if not nb_samples:
    nb_samples = 4

rule spliceLauncher_sashimi_plot:
    input:
        bam = os.path.abspath(path_bam + "spliceLauncher_STAR/" + name_genome + "/{reads}.sorted.bam"),
        event_file = path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.{junction}.filter" + extension,

    output:
        directory = directory(path_results + "/spliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/sashimi_plot/{junction}/")

    params:
        list_bam = expand("{reads}", reads = all_samples),
        python = python,
        ggsashimi = ggsashimi,
        extend_bp = extend_bp,
        MinThresholdNbReads = MinThresholdNbReads,
        nb_samples = nb_samples,
        gtf = gtf,
        color = color,

    shell:
        "{params.python} scripts/generate_sashimi_plot.py "
        "-directory {output.directory} "
        "-ggsashimi {params.ggsashimi} "
        "-bam {input.bam} "
        "-event_file {input.event_file} "
        "-extend_bp {params.extend_bp} "
        "-MinThresholdNbReads {params.MinThresholdNbReads} "
        '-list_bam "{params.list_bam}" '
        "-nb_samples {params.nb_samples} "
        "-gtf {params.gtf} "
        "-color {params.color}"





    #     parser.add_argument('-ggsashimi', dest='ggsashimi', help='full path to executable ggsashimi')
    # parser.add_argument('-bam', dest='bam', help='full path to interest bam')
    # parser.add_argument('-event_file', dest='event_file', help='full path to event interest')
    # parser.add_argument('-extend_bp', dest='extend_bp', help='number of bp add to extremity. Default = 50', type=int, default=50)
    # parser.add_argument('-MinThresholdNbReads', dest='MinThresholdNbReads', help='Minimum reads in junctions for plot. Default, half of the reads of the sample of interest junction', type=int, default=-1)
    # parser.add_argument('-list_bam', dest='list_bam', help='list of bam')
    # parser.add_argument('-nb_samples', dest='nb_samples', help='number of random sample in plot. Default = 4', type=int, default=4)
    # parser.add_argument('-gtf', dest='gtf', help='full path to gtf of reference', default=-1)
    # parser.add_argument('-color', dest='color', help='palette of color for plot', default=-1)
    # args = parser.parse_args()