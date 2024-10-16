#Quality controle of raw reads
rule fastqc_raw:
    input:
        read = os.path.abspath(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_process/{reads}.fastq.gz")

    output:
        html = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_raw/{reads}_fastqc.html",
        zip = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_raw/{reads}_fastqc.zip"


    params:
        fastqc = config["DEPENDANCES"]["FASTQC"],
        directory = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_raw/"

    threads:
        config["PARAMS"]["FASTQC"]["THREADS"]
    
    log:
        out = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/fastqc_raw/{reads}.stdout.log",
        err = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/fastqc_raw/{reads}.stderr.log"

    run:
        create_directory_if_not_exists(params["directory"])
        shell("{params.fastqc} "
        "-o {params.directory} "
        "-t {threads} "
        "{input.read} > {log.out} 2> {log.err}")



#Quality controle of trimmed reads
rule fastqc_trimmed:
    input:
        read = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_trimmed/{reads}.trimmed.fastq.gz"

    output:
        html = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_trimmed/{reads}.trimmed_fastqc.html",
        zip = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_trimmed/{reads}.trimmed_fastqc.zip"


    params:
        fastqc = config["DEPENDANCES"]["FASTQC"],
        directory = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_trimmed/"

    threads:
        config["PARAMS"]["FASTQC"]["THREADS"]
    
    log:
        out = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/fastqc_trimmed/{reads}.stdout.log",
        err = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/fastqc_trimmed/{reads}.stderr.log"

    run:
        create_directory_if_not_exists(params["directory"])
        shell("{params.fastqc} "
        "-o {params.directory} "
        "-t {threads} "
        "{input.read} > {log.out} 2> {log.err}")


#Rapport of fastq
rule multiqc_fastq_raw:
    input:
        html1 = expand(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_raw/{reads}_1_fastqc.html",reads = all_samples),
        zip1 = expand(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_raw/{reads}_1_fastqc.zip", reads = all_samples),
        html2 = expand(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_raw/{reads}_2_fastqc.html",reads = all_samples),
        zip2 = expand(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_raw/{reads}_2_fastqc.zip", reads = all_samples)

    output:
        directory_data = directory(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/multiqc/" + config["PARAMS"]["GENERAL"]["PREFIX"] +"_fastq_raw_" + unique_id + "_data/"),
        html =config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/multiqc/" + config["PARAMS"]["GENERAL"]["PREFIX"]  + "_fastq_raw_" + unique_id + ".html"
        
    params:
        name = config["PARAMS"]["GENERAL"]["PREFIX"] + "_fastq_raw_" + unique_id,
        multiqc = config["DEPENDANCES"]["MULTIQC"],
        path = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/multiqc/"

    threads:
        config["PARAMS"]["MULTIQC"]["THREADS"]
        
    shell:
        "{params.multiqc} "
        "{input.html1} {input.zip1} {input.html2} {input.zip2} "
        "--filename {params.name} "
        "-o {params.path}"


#Rapport of fastq
rule multiqc_fastq_trimmed:
    input: 
        html1 = expand(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_trimmed/{reads}_1.trimmed_fastqc.html",reads = all_samples),
        zip1 = expand(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_trimmed/{reads}_1.trimmed_fastqc.zip", reads = all_samples),
        html2 = expand(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_trimmed/{reads}_2.trimmed_fastqc.html",reads = all_samples),
        zip2 = expand(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc_trimmed/{reads}_2.trimmed_fastqc.zip", reads = all_samples),
        html = expand(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastp/{reads}_fastp.html", reads = all_samples),
        json = expand(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastp/{reads}_fastp.json", reads = all_samples)

    output:
        directory_data = directory(config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/multiqc/" + config["PARAMS"]["GENERAL"]["PREFIX"] +"_fastq_trimmed_" + unique_id + "_data/"),
        html = config["DATA_INPUTS"]["WORKING_DIRECTORY"]  + "/2-processed_data/quality_control/multiqc/" + config["PARAMS"]["GENERAL"]["PREFIX"]  + "_fastq_trimmed_" + unique_id + ".html"

    params:
        name = config["PARAMS"]["GENERAL"]["PREFIX"] + "_fastq_trimmed_" + unique_id,
        multiqc = config["DEPENDANCES"]["MULTIQC"],
        path = config["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/multiqc/"

    threads:
        config["PARAMS"]["MULTIQC"]["THREADS"]
        
    shell:
        "{params.multiqc} "
        "{input.html1} {input.zip1} {input.html2} {input.zip2} {input.html} {input.json} "
        "--filename {params.name} "
        "-o {params.path}"