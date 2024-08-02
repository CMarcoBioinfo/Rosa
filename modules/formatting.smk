# Compresses a read file
rule compress_fastq:
    input:
        read = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/{reads}"

    output:
        read_final = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/{reads}.gz"
    
    params:
        pigz = config["DEPENDANCES"]["PIGZ"]

    threads: 
        config["PARAMS"]["COMPRESS"]["THREADS"]

    priority: 3

    shell:
        "{params.pigz} "
        "--best "
        "--processes {threads} {input.read}"


#Renames a read file
rule rename_reads:
    input:
        read = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/reads/{read}.fq.gz"
    
    output:
        read_final = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/reads/{read}.fastq.gz"
    
    priority: 2

    shell:
        "mv "
        "{input.read} "
        "{output.read_final}"


#Deinterleaves a fastq file into two separate fastq file
rule deinterleave_fastq:
    input:
        read = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/reads/{reads}.fastq.gz"
    
    output:
        read1 = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/samples/reads1/{reads}_1.fastq.gz",
        read2 = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/samples/reads2/{reads}_2.fastq.gz"

    params:
        pigz = config["DEPENDANCES"]["PIGZ"]

    threads:
        config["PARAMS"]["COMPRESS"]["THREADS"]

    priority: 1

    shell:
        '{params.pigz} --best --processes {threads} -dc {input.read} | '
        'paste -------- | '
        'tee > '
        '(cut -f 1-4 | tr "\\t" "\\n" | pigz --best --processes {threads} > {output.read1}) | '
        'cut -f 5-8 | tr "\\t" "\\n" | pigz --best --processes {threads} > {output.read2}'


#Sort BAM file
rule samtools_sort:
    input:
        bam = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/bam/{reads}.bam"
    
    output:
        bam_sort = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/bam/{reads}_sorted.bam"
    
    params:
        samtools = config["DEPENDANCES"]["SAMTOOLS"]
    
    threads: 
        config["PARAMS"]["SAMTOOLS"]["THREADS"]

    shell:
        "{params.samtools} sort "
        "-n {input.bam} "
        "-@ {threads} "
        "-o {output.bam_sort}"


#Converts a sorted BAM file into two separate FASTQ files.
rule samtools_fastq:
    input: 
        bam_sort = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/bam/{reads}_sorted.bam"
    
    output:
        read1 = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/samples/reads1/{reads}_1.fastq.gz",
        read2 = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/samples/reads2/{reads}_2.fastq.gz"
    
    params:
        samtools = config["DEPENDANCES"]["SAMTOOLS"]
    
    threads: 
        config["PARAMS"]["SAMTOOLS"]["THREADS"]
    
    shell:
        "{params.samtools} fastq "
        "@ {threads} "
        "{input.bam_sort} "
        "-1 {output.read1} "
        "-2 {output.read2} "
        "-0 /dev/null -s /dev/null -n "
        "&& rm {input.bam_sort}"