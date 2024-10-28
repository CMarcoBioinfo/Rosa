#hisat2-build -p 3 -f GRCh37.p13.chr.fa  GRCh37.p13

# hisat2 -p 3 --dta -x ../genome/GRCh37.p13 -1 ../../../data/01-ROSALIND-fastq/202304-1409241167-SLA_1.fastq.gz -2 ../../../data/01-ROSALIND-fastq/202304-1409241167-SLA_2.fastq.gz | samtools view -bS - | samtools sort -T tmp/202304-1409241167-SLA/ -o output/202304-1409241167-SLA/202304-1409241167-SLA.sorted.bam -

#rule compress_fastq:
#    input:



#rule rename_reads:
#    input:


#rule deinterleave_fastq:
#    input:



#Index reference
rule hisat_index:
    input:
        reference = os.path.abspath(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"])

    output:
        directory = directory(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0]),
        reference = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"]
        
    params:
        hisat2 = config["DEPENDANCES"]["MAPPING"]["HISAT2"],
        name_genome = "/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0]
        
    threads:
        config["PARAMS"]["HISAT2"]["THREADS"]
    
    log:
        out = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/hisat_index/"+ config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] +".stdout.log",
        err = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/hisat_index/"+ config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] +".stderr.log"


    run:
        create_directory_if_not_exists(output["directory"])
        shell("ln -sfn {input.reference} {output.reference} && "
        "{params.hisat2}-build "
        "-p {threads} "
        "{output.reference} {output.directory}{params.name_genome} > {log.out} 2> {log.err}")

    
    

#Align reads to reference
rule hisat_align:
    input:
        directory = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0],
        reference = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0] + "/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"],
        # reads1 = lambda wcs: reads(wcs, "1"),
        # reads2 =  lambda wcs: reads(wcs, "2")
        trimmed_read1 = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_trimmed/{reads}_1.trimmed.fastq.gz",
        trimmed_read2 = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/2-processed_data/samples_trimmed/{reads}_2.trimmed.fastq.gz"

    output:
        bam = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam",
        summary = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.summary.txt"


    params:
        name_genome = "/" + config["GENERAL"]["DATA_INPUTS"]["GENOME"].rsplit(".",1)[0],
        hisat2 = config["DEPENDANCES"]["MAPPING"]["HISAT2"],
        samtools = config["DEPENDANCES"]["MAPPING"]["SAMTOOLS"]
    
    log:
        out = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/hisat2_align/{reads}.stdout.log",
        err = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/hisat2_align/{reads}.stderr.log"

    threads:
        config["PARAMS"]["HISAT2"]["THREADS"]

    shell:
        "{params.hisat2} "
        "-p {threads} "
        "--summary-file {output.summary} "
        "-x {input.directory}{params.name_genome} "
        "-1 {input.trimmed_read1} "
        "-2 {input.trimmed_read2} | "
        "{params.samtools} sort "
        "-@ {threads} "
        "-o {output.bam} - > {log.out} 2> {log.err}"


#Index BAM
rule samtools_index:
    input:
        bam = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam"

    output:
        bai = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam.bai"
    
    params:
        wd = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"],
        sample = "{reads}",
        samtools = config["DEPENDANCES"]["MAPPING"]["SAMTOOLS"]
    
    log:
        out = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/samtools_index/{reads}.stdout.log",
        err = config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/log/samtools_index/{reads}.stderr.log"
    
    priority: 1 
    
    shell:
        "{params.samtools} index "
        "{input.bam} > {log.out} 2> {log.err}"