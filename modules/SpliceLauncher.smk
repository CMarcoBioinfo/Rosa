rule SpliceLauncher_create_bed:
    input:
        bam = os.path.abspath(path_bam + name_genome + "/mapping/{reads}.sorted.bam"),
        csi = path_bam + name_genome + "/mapping/{reads}.sorted.bam.csi"

    output:
        bed = path_bam  + name_genome + "/SpliceLauncher/{reads}_juncs.bed"

    params:
        samtools = samtools,
        bedtools = bedtools,
        SpliceLauncher = SpliceLauncher,
        perl = perl

    shell:
        "{params.samtools} view "
        "-b -f 0x40 -F 1024 {input.bam} | "
        "{params.bedtools} bamtobed -bed12 -i stdin | "
        "awk '{{if($10>1){{print $0}}}}' | "
        "{params.perl} {params.SpliceLauncher}/scripts/bedBlocks2IntronsCoords.pl y - | "
        "awk '{{if($5==255){{print $0}}}}' > {output.bed} && "
        "{params.samtools} view "
        "-b -f 0x80 -F 1024 {input.bam} | "
        "{params.bedtools} bamtobed -bed12 -i stdin | "
        "awk '{{if($10>1){{print $0}}}}' | "
        "{params.perl} {params.SpliceLauncher}/scripts/bedBlocks2IntronsCoords.pl n - | "
        "awk '{{if($5==255){{print $0}}}}' >> {output.bed}"


rule SpliceLauncher_count_Junctions:
    input:
        bed_ref = working_directory + "/2-processed_data/reference/" + name_genome + "/BEDannotation.bed",
        bed = path_bam  + name_genome + "/SpliceLauncher/{reads}_juncs.bed"

    output:
        counts = path_bam  + name_genome + "/SpliceLauncher/getClosestExons/{reads}.count"

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


rule SpliceLauncher_merge_count:
    input:
        counts = expand(path_bam  + name_genome + "/SpliceLauncher/getClosestExons/{reads}.count",reads= all_samples)

    output:
        mergeFile = path_results + "/SpliceLauncher/mergeFile/" + prefix + "_" + unique_id + ".txt"
    
    params:
        SpliceLauncher = SpliceLauncher,
        perl = perl,
        directory = directory(path_results + "/SpliceLauncher/mergeFiles/"),
        genes_of_interest = genes_of_interest

    shell:
        "{params.perl} {params.SpliceLauncher}/scripts/joinJuncFiles.pl -c {input.counts} {params.genes_of_interest} > {output.mergeFile}"


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
count_report = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_report_" + date + ".txt",
listOutputSpliceLauncher.append(count_report)
outputSpliceLauncher = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + "_outputSpliceLauncher" + extension,
listOutputSpliceLauncher.append(outputSpliceLauncher)

bedOut = config["SPLICELAUNCHER"]["ANALYSE"].get("BED_OUT")
if bedOut:
    bedOut = "--bedOut "
    bed = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + ".bed"
    listOutputSpliceLauncher.append(bed)

else:
    bedOut = ""

Graphics = config["SPLICELAUNCHER"]["ANALYSE"].get("GRAPHICS")
if Graphics:
    Graphics = "--Graphics "
    pdf = expand(path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.pdf",reads= all_samples)
    pdf_significant_genes = expand(path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.significant_genes.pdf",reads= all_samples)
    pdf_significant_junctions = expand(path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.significant_junctions.pdf",reads= all_samples)
    listOutputSpliceLauncher.append(pdf)
    listOutputSpliceLauncher.append(pdf_significant_genes)
    listOutputSpliceLauncher.append(pdf_significant_junctions)




rule SpliceLauncher_Analyse:
    input:
        mergeFile = path_results + "/SpliceLauncher/mergeFile/" + prefix + "_" + unique_id + ".txt",
        annot = working_directory + "/2-processed_data/reference/" + name_genome + "/SpliceLauncherAnnot.txt"

    output:
        listOutputSpliceLauncher

    params:
        SpliceLauncher = SpliceLauncher,
        Rscript = Rscript,
        nbIntervals = nbIntervals,
        threshold = threshold,
        min_cov = min_cov,
        transcriptList = transcriptList,
        removeOther = removeOther,
        txt = txt,
        bedOut = bedOut,
        Graphics = Graphics,
        directory = path_results + "/SpliceLauncher/",
        sampleNames = sampleNames


    shell:
        "{params.Rscript} {params.SpliceLauncher}/scripts/SpliceLauncherAnalyse.r "
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
filterFilesSampleSignificant = expand(path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.significant_junctions" + extension,reads= all_samples)

filterFileSignificant = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + "_outputSpliceLauncher.significant_junctions" + extension,
filterFileUnique = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + "_outputSpliceLauncher.unique_junctions" + extension,
outputFilterAnalyse.append(filterFileSignificant)
outputFilterAnalyse.append(filterFileUnique)
outputFilterAnalyse.append(filterFilesSampleSignificant)

filterFilesSampleSignificantFilter = expand(path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.significant_junctions.filter" + extension,reads= all_samples)
filterFilesSampleUniqueFilter = expand(path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.unique_junctions.filter" + extension,reads= all_samples)
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


rule SpliceLauncher_filter_analyse:
    input:
        outputSpliceLauncher = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/" +  prefix + "_" + unique_id + "_outputSpliceLauncher" + extension,

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
        directory = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/"

    shell:
        "{params.Rscript} scripts/SpliceLauncher_filter_analyse.r "
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

rule SpliceLauncher_sashimi_plot:
    input:
        bam = os.path.abspath(path_bam + name_genome + "/mapping/{reads}.sorted.bam"),
        event_file = path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/{reads}.{junction}.filter" + extension,

    output:
        directory = directory(path_results + "/SpliceLauncher/" + prefix + "_" + unique_id + "_results/samples_results/{reads}/sashimi_plot/{junction}/")

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