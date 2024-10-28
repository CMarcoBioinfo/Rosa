# Compresses a read file
rule compress_fastq:
    input:
        read = lambda wildcards: get_inputs(wildcards,0),
        preprocess = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/.tmp/samples/{id}.pre"


    output:
        process = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/.tmp/samples/{id}.pro"

    
    params:
        directory = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_raw/",
        id = lambda wildcards: wildcards.id,
        pigz = config["DEPENDANCES"]["FORMATING"]["PIGZ"]

    threads: 
        config["PARAMS"]["COMPRESS"]["THREADS"]
    
    priority: 5

    run:
        if (str(input.read).endswith('.gz') ):
            shell("touch {output.process}")
        else : 
            shell("{params.pigz} "
            "--best "
            "--verbose "
            "--processes {threads} {input.read} && touch {output.process}")



# #Deinterleaves a fastq file into two separate fastq file
# rule deinterleave_fastq:
#     input:
#         read = lambda wildcards: get_inputs(wildcards,3),
#         process = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/.tmp/samples/reads/{id}.pro"
    
#     output:
#         read1 = os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_process/{id}_1.fastq.gz"),
#         read2 = os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_process/{id}_2.fastq.gz")

#     params:
#         id = lambda wildcards: wildcards.id,
#         pigz = config["DEPENDANCES"]["FORMATING"]["PIGZ"],
#         directory = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_process/"

#     priority: 6
    
#     threads:
#         config["PARAMS"]["COMPRESS"]["THREADS"]

#     run:
#         create_directory_if_not_exists(params["directory"])
#         shell('{params.pigz} --best --processes {threads} -dc {input.read} | '
#         'paste -------- | '
#         'tee > '
#         '(cut -f 1-4 | tr "\\t" "\\n" | pigz --best --processes {threads} > {output.read1}) | '
#         'cut -f 5-8 | tr "\\t" "\\n" | pigz --best --processes {threads} > {output.read2}')


# #Sort BAM file
# rule samtools_sort:
#     input:
#         bam = lambda wildcards: get_inputs(wildcards,4),
#         preprocess = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/.tmp/bam/{id}.pre"

    
#     output:
#         bam_sort = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_raw/bam/{id}_sorted.bam",
#         process = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/.tmp/bam/{id}.pro"

    
#     params:
#         id = lambda wildcards: wildcards.id,
#         samtools = config["DEPENDANCES"]["SAMTOOLS"],
#         directory = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_raw/bam"

#     threads: 
#         config["PARAMS"]["SAMTOOLS"]["THREADS"],

#     log:
#         out = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/samtools_sort/{id}.stout.log",
#         err = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/samtools_sort/{id}.stderr.log"
    

#     run:
#         create_directory_if_not_exists(params["directory"])
#         shell("{params.samtools} sort "
#         "-n {input.bam} "
#         "-@ {threads} "
#         "-o {output.bam_sort} > {log.out} 2> {log.err} && touch {output.process}")


# #Converts a sorted BAM file into two separate FASTQ files.
# rule samtools_fastq:
#     input: 
#         bam_sort = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_raw/bam/{id}_sorted.bam",
#         process = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/.tmp/bam/{id}.pro"

    
#     output:
#         read1 = os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_process/{id}_1.fastq.gz"),
#         read2 = os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_process/{id}_2.fastq.gz")
    
#     params:
#         samtools = config["DEPENDANCES"]["SAMTOOLS"],
#         directory = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_process/"
    
#     threads: 
#         config["PARAMS"]["SAMTOOLS"]["THREADS"]

#     priority: 1

#     run:
#         create_directory_if_not_exists(params["directory"])
#         shell("{params.samtools} fastq "
#         "-@ {threads} "
#         "-1 {output.read1} "
#         "-2 {output.read2} "
#         "-0 /dev/null -s /dev/null -n {input.bam_sort} "
#         "&& rm {input.bam_sort}")

rule processed_fastq:
    input:
        read1 = lambda wildcards: get_inputs(wildcards,1),
        read2 =  lambda wildcards: get_inputs(wildcards,2),
        process1 = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/.tmp/samples/reads1/{id}.pro",
        process2 = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/.tmp/samples/reads2/{id}.pro"

    output:
        read1 = os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_process/{id}_1.fastq.gz"),
        read2 = os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_process/{id}_2.fastq.gz")

    
    params:
        id = lambda wildcards: wildcards.id,
        directory = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_process/"

    priority: 5

    run:
        create_directory_if_not_exists(params["directory"])
        if (is_gz(input["read1"]) == False):
            input["read1"] = input["read1"] + ".gz"
        if (is_gz(input["read2"]) == False):
            input["read2"] = input["read2"] + "gz"
        shell("ln -sfn "
        "{input.read1} {output.read1} && "
        "ln -sfn "
        "{input.read2} {output.read2}")


#Trimming adaptator from raw reads
rule fastp_trimming:
    input:
        read1 = os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_process/{reads}_1.fastq.gz"),
        read2 = os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_process/{reads}_2.fastq.gz")

    output:
        trimmed_read1 = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_trimmed/{reads}_1.trimmed.fastq.gz",
        trimmed_read2 = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_trimmed/{reads}_2.trimmed.fastq.gz",
        html = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastp/{reads}_fastp.html",
        json = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastp/{reads}_fastp.json"

    params:
        fastp = config["DEPENDANCES"]["FORMATING"]["FASTP"],
        length = config["PARAMS"]["FASTP"]["LENGTH"],
        directory_trimmed = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_trimmed/",
        directory_fastp = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastp/"


    threads:
        config["PARAMS"]["FASTP"]["THREADS"]

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


rule gtf2bed:
    input:
        gtf = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/annotation/" + config["GENERAL"]["DATA_INPUTS"]["ANNOTATION"]

    output:
        bed = (config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/annotation/" + config["GENERAL"]["DATA_INPUTS"]["ANNOTATION"]).rsplit(".",1)[0] + ".bed",
        bed3 = (config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/annotation/" + config["GENERAL"]["DATA_INPUTS"]["ANNOTATION"]).rsplit(".",1)[0] + ".bed3"

    params:
        gtf2bed = config["DEPENDANCES"]["FORMATING"]["PIGZ"]

    shell:
        "gtf2bed " 
        "< {input.gtf} "
        "> {output.bed} "