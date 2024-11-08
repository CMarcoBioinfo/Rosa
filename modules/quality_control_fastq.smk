#Quality controle of raw reads
rule fastqc_raw:
    input:
        read = os.path.abspath(working_directory + "/1-raw_data/fastq/{reads}.fastq.gz")

    output:
        html = path_qc + "fastqc_raw/{reads}_fastqc.html",
        zip = path_qc + "fastqc_raw/{reads}_fastqc.zip"


    params:
        fastqc = fastqc,
        directory = path_qc + "fastqc_raw/"

    threads:
        config["PARAMS"]["FASTQC"]["THREADS"]
    
    log:
        out = working_directory  + "/log/fastqc_raw/{reads}.out",
        err = working_directory  + "/log/fastqc_raw/{reads}.err"

    run:
        create_directory_if_not_exists(params["directory"])
        shell("{params.fastqc} "
        "-o {params.directory} "
        "-t {threads} "
        "{input.read} > {log.out} 2> {log.err}")


#Rapport of fastq
rule multiqc_fastq_raw:
    input:
        html1 = expand(path_qc + "fastqc_raw/{reads}.1_fastqc.html",reads = all_samples),
        zip1 = expand(path_qc + "fastqc_raw/{reads}.1_fastqc.zip", reads = all_samples),
        html2 = expand(path_qc + "fastqc_raw/{reads}.2_fastqc.html",reads = all_samples),
        zip2 = expand(path_qc + "fastqc_raw/{reads}.2_fastqc.zip", reads = all_samples)

    output:
        directory_data = directory(path_qc + "multiqc/fastq_raw/" + prefix +"_" + unique_id + "_data/"),
        html = path_qc + "multiqc/fastq_raw/" + prefix + "_" + unique_id + ".html"
        
    params:
        name = prefix + "_" + unique_id,
        multiqc = multiqc,
        path = path_qc + "multiqc/fastq_raw/"

    threads:
        config["PARAMS"]["MULTIQC"]["THREADS"]
        
    shell:
        "{params.multiqc} "
        "{input.html1} {input.zip1} {input.html2} {input.zip2} "
        "--filename {params.name} "
        "-o {params.path}"

#Quality controle of trimmed reads
rule fastqc_trimmed:
    input:
        trimmed_read = os.path.abspath(path_fastq + "{reads}.fastq.gz")

    output:
        html = path_qc + "fastqc_trimmed/{reads}_fastqc.html",
        zip = path_qc + "fastqc_trimmed/{reads}_fastqc.zip"


    params:
        fastqc = fastqc,
        directory = path_qc + "fastqc_trimmed/"

    threads:
        config["PARAMS"]["FASTQC"]["THREADS"]
    
    log:
        out = working_directory  + "/log/fastqc_trimmed/{reads}.out",
        err = working_directory  + "/log/fastqc_trimmed/{reads}.err"

    run:
        create_directory_if_not_exists(params["directory"])
        shell("{params.fastqc} "
        "-o {params.directory} "
        "-t {threads} "
        "{input.trimmed_read} > {log.out} 2> {log.err}")


#Rapport of fastq
rule multiqc_fastq_trimmed:
    input: 
        html1 = expand(path_qc + "fastqc_trimmed/{reads}.1_fastqc.html",reads = all_samples),
        zip1 = expand(path_qc + "fastqc_trimmed/{reads}.1_fastqc.zip", reads = all_samples),
        html2 = expand(path_qc + "fastqc_trimmed/{reads}.2_fastqc.html",reads = all_samples),
        zip2 = expand(path_qc + "fastqc_trimmed/{reads}.2_fastqc.zip", reads = all_samples),
        html = expand(path_qc + "fastp_trimming/{reads}.html", reads = all_samples),
        json = expand(path_qc + "fastp_trimming/{reads}.json", reads = all_samples)

    output:
        directory_data = directory(path_qc + "multiqc/fastq_trimmed/" + prefix +"_" + unique_id + "_data/"),
        html = path_qc + "multiqc/fastq_trimmed/" + prefix +"_" + unique_id + ".html"

    params:
        name = prefix +"_" + unique_id,
        multiqc = multiqc,
        path =  path_qc + "multiqc/fastq_trimmed/"

    threads:
        config["PARAMS"]["MULTIQC"]["THREADS"]
        
    shell:
        "{params.multiqc} "
        "{input.html1} {input.zip1} {input.html2} {input.zip2} {input.html} {input.json} "
        "--filename {params.name} "
        "-o {params.path}"