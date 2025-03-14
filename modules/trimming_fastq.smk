#Trimming adaptator from raw reads

rule fastp_trimming:
    input:
        read1 = os.path.abspath(working_directory + "/1-raw_data/fastq/{reads}.1.fastq.gz"),
        read2 = os.path.abspath(working_directory + "/1-raw_data/fastq/{reads}.2.fastq.gz")

    output:
        trimmed_read1 = os.path.abspath(path_fastq + "{reads}.1.fastq.gz"),
        trimmed_read2 = os.path.abspath(path_fastq + "{reads}.2.fastq.gz"),
        html = path_qc + "fastp_trimming/{reads}.html",
        json = path_qc + "fastp_trimming/{reads}.json"

    params:
        fastp = fastp,
        length = length_fastp,
        directory_trimmed = path_fastq,
        directory_fastp = path_qc + "fastp_trimming/"


    threads:
        config["GENERAL"]["THREADS"]

    run:
        create_directory_if_not_exists(params["directory_trimmed"])
        create_directory_if_not_exists(params["directory_fastp"])
        shell("{params.fastp} "
        "--thread {threads} "
        "--length_required {params.length} "
        "--detect_adapter_for_pe "
        "--in1 {input.read1} "
        "--in2 {input.read2} "
        "--out1 {output.trimmed_read1} "
        "--out2 {output.trimmed_read2} "
        "--html {output.html} "
        "--json {output.json} ")