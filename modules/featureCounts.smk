from scripts import merge_counts as mc

rule featureCounts:
    input:
        bai = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam.bai",
        annotation = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/annotation/" + config["DATA_INPUT"]["ANNOTATION"]

    output:
        summaryFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/{reads}.counts.summary",
        countsFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/{reads}.counts",


    params:
        directory = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/",
        bam = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam",
        featureCounts = config["DEPENDANCES"]["FEATURECOUNTS"],
        reads = "{reads}",
        metadata = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/sample_names.txt"



    threads:
        config["PARAMS"]["FEATURECOUNTS"]["THREADS"]

    run:
        create_directory_if_not_exists(params["directory"])
        shell("{params.featureCounts} -p "
        "-T {threads} "
        "-a {input.annotation} "
        "-o {params.directory}{params.reads}.counts {params.bam}")
        add_file(params["reads"], params["bam"], output.countsFile)


rule merge_counts:
    input:
        counts = expand(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/{reads}.counts", reads = all_samples),

    output:
        mergeFile = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/merge.counts"
    
    params:
        metadata = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/2-Counts/sample_names.txt"

    threads:
        config["PARAMS"]["MERGING"]["THREADS"]


    run:
        mc.merging_files(params["metadata"],output.mergeFile, threads)


