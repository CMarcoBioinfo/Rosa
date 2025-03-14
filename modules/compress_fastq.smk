# Compresses a read file
rule compress_fastq:
    input:
        read = lambda wildcards: get_inputs(wildcards,0),
        preprocess = working_directory + "/.tmp/samples/{id}.pre"
        
    output:
        process = working_directory + "/.tmp/samples/{id}.pro"

    
    params:
        id = lambda wildcards: wildcards.id,
        pigz = pigz

    threads: 
        config["GENERAL"]["THREADS"]
    
    priority: 5

    run:
        if (str(input.read).endswith('.gz') ):
            if (not os.path.exists(output.process)):
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
#         process = working_directory + "/.tmp/samples/reads/{id}.pro"
    
#     output:
#         read1 = os.path.abspath(working_directory + "/2-processed_data/samples_process/{id}_1.fastq.gz"),
#         read2 = os.path.abspath(working_directory + "/2-processed_data/samples_process/{id}_2.fastq.gz")

#     params:
#         id = lambda wildcards: wildcards.id,
#         pigz = pigz,
#         directory = working_directory + "/2-processed_data/samples_process/"

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
#         preprocess = working_directory + "/.tmp/bam/{id}.pre"

    
#     output:
#         bam_sort = working_directory + "/2-processed_data/samples_raw/bam/{id}_sorted.bam",
#         process = working_directory + "/.tmp/bam/{id}.pro"

    
#     params:
#         id = lambda wildcards: wildcards.id,
#         samtools = config["DEPENDANCES"]["SAMTOOLS"],
#         directory = working_directory + "/2-processed_data/samples_raw/bam"

#     threads: 
#         config["PARAMS"]["SAMTOOLS"]["THREADS"],

#     log:
#         out = working_directory + "/log/samtools_sort/{id}.stout.log",
#         err = working_directory + "/log/samtools_sort/{id}.stderr.log"
    

#     run:
#         create_directory_if_not_exists(params["directory"])
#         shell("{params.samtools} sort "
#         "-n {input.bam} "
#         "-@ {threads} "
#         "-o {output.bam_sort} > {log.out} 2> {log.err} && touch {output.process}")


# #Converts a sorted BAM file into two separate FASTQ files.
# rule samtools_fastq:
#     input: 
#         bam_sort = working_directory + "/2-processed_data/samples_raw/bam/{id}_sorted.bam",
#         process = working_directory + "/.tmp/bam/{id}.pro"

    
#     output:
#         read1 = os.path.abspath(working_directory + "/2-processed_data/samples_process/{id}_1.fastq.gz"),
#         read2 = os.path.abspath(working_directory + "/2-processed_data/samples_process/{id}_2.fastq.gz")
    
#     params:
#         samtools = config["DEPENDANCES"]["SAMTOOLS"],
#         directory = working_directory + "/2-processed_data/samples_process/"
    
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
        process1 = working_directory + "/.tmp/samples/reads1/{id}.pro",
        process2 = working_directory + "/.tmp/samples/reads2/{id}.pro"

    output:
        read1 = os.path.abspath(working_directory + "/1-raw_data/fastq/{id}.1.fastq.gz"),
        read2 = os.path.abspath(working_directory + "/1-raw_data/fastq/{id}.2.fastq.gz")

    
    params:
        id = lambda wildcards: wildcards.id,
        directory = working_directory + "/1-raw_data/fastq/"

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




# rule gtf2bed:
#     input:
#         gtf = working_directory + "/1-raw_data/annotation/" + config["GENERAL"]["ANNOTATION"]

#     output:
#         bed = (working_directory + "/1-raw_data/annotation/" + config["GENERAL"]["ANNOTATION"]).rsplit(".",1)[0] + ".bed",
#         bed3 = (working_directory + "/1-raw_data/annotation/" + config["GENERAL"]["ANNOTATION"]).rsplit(".",1)[0] + ".bed3"

#     params:
#         gtf2bed = config["DEPENDANCES"]["FORMATING"]["PIGZ"]

#     shell:
#         "gtf2bed " 
#         "< {input.gtf} "
#         "> {output.bed} "