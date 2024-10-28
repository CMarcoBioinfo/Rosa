from scripts import merge_counts as mc

rule featureCounts:
    input:
        bai = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam.bai",
        annotation = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/annotation/" + config["GENERAL"]["DATA_INPUTS"]["ANNOTATION"]

    output:
        summaryFile = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/2-Counts/{reads}.counts.summary",
        countsFile = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/2-Counts/{reads}.counts",


    params:
        directory = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/2-Counts/",
        bam = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam",
        featureCounts = config["DEPENDANCES"]["ANALYSES"]["FEATURECOUNTS"],
        reads = "{reads}",
        metadata = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/2-Counts/sample_names.txt"



    threads:
        config["PARAMS"]["FEATURECOUNTS"]["THREADS"]

    run:
        create_directory_if_not_exists(params["directory"])
        shell("{params.featureCounts} -p --countReadPairs "
        "-T {threads} "
        "-a {input.annotation} "
        "-o {params.directory}{params.reads}.counts {params.bam}")
        add_file(params["reads"], params["bam"], output.countsFile)


rule merge_counts:
    input:
        counts = expand(config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/2-Counts/{reads}.counts", reads = all_samples),

    output:
        mergeFile = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/2-Counts/merge.counts"
    
    params:
        metadata = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/2-Counts/sample_names.txt"

    threads:
        config["PARAMS"]["MERGING"]["THREADS"]

    run:
        mc.merging_files(params["metadata"],output.mergeFile, threads)


