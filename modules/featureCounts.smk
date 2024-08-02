# 1997  featureCounts -p -T 2 -a ../../../../06-data/05-GTF/GRCh37.P13/GRCh37.p13.SLA.chr.gtf -o test.txt 202304-1409241167-SLA.sorted.bam
#Include script
from scripts import merge_counts as mc

rule featureCounts:
    input:
        bai = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam.bai",
        annotation = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/annotation/" + config["DATA_INPUT"]["ANNOTATION"]

    output:
        summaryFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/{reads}.counts.summary",
        countsFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/{reads}.counts"

    params:
        directory = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/",
        bam = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam",
        featureCounts = config["DEPENDANCES"]["FEATURECOUNTS"],
        reads = "{reads}"

    threads:
        config["PARAMS"]["FEATURECOUNTS"]["THREADS"]

    run:
        create_directory_if_not_exists(params["directory"])
        shell("{params.featureCounts} -p "
        "-T {threads} "
        "-a {input.annotation} "
        "-o {params.directory}{params.reads}.counts {params.bam}")


rule merge_file:
    input:
        countsFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/{reads}.counts"
    
    output:
        controlFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/.Control/{reads}.txt"
    
    params:
        mergeFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/merge.counts",
        directory = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/.Control/",


    threads:
        config["PARAMS"]["MERGING"]["THREADS"]
    
    run:
        create_directory_if_not_exists(params["directory"])
        mc.merging_files(params.mergeFile, input.countsFile, output.controlFile, threads)