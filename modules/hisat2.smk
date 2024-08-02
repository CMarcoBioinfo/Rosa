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
        reference = os.path.abspath(config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/reference/" + config["DATA_INPUT"]["GENOME"])

    output:
        directory = directory(config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["DATA_INPUT"]["GENOME"].rsplit(".",1)[0]),
        reference = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["DATA_INPUT"]["GENOME"].rsplit(".",1)[0] + "/" + config["DATA_INPUT"]["GENOME"]
        
    params:
        hisat2 = config["DEPENDANCES"]["HISAT2"],
        name_genome = "/" + config["DATA_INPUT"]["GENOME"].rsplit(".",1)[0]
        

    threads:
        config["PARAMS"]["HISAT2"]["THREADS"]

    run:
        create_directory_if_not_exists(output["directory"])
        shell("ln -sfn {input.reference} {output.reference} && "
        "{params.hisat2}-build "
        "-p {threads} "
        "{output.reference} {output.directory}{params.name_genome}")

    
    

#Align reads to reference
rule hisat_align:
    input:
        directory = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["DATA_INPUT"]["GENOME"].rsplit(".",1)[0],
        reference = config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/reference/" + config["DATA_INPUT"]["GENOME"].rsplit(".",1)[0] + "/" + config["DATA_INPUT"]["GENOME"],
        reads1 = lambda wcs: reads(wcs, "1"),
        reads2 =  lambda wcs: reads(wcs, "2")

    output:
        bam = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam"

    params:
        name_genome = "/" + config["DATA_INPUT"]["GENOME"].rsplit(".",1)[0],
        hisat2 = config["DEPENDANCES"]["HISAT2"],
        samtools = config["DEPENDANCES"]["SAMTOOLS"]

    threads:
        config["PARAMS"]["HISAT2"]["THREADS"]

    shell:
        "{params.hisat2} "
        "-p {threads} "
        "-x {input.directory}{params.name_genome} "
        "-1 {input.reads1} "
        "-2 {input.reads2} | "
        "{params.samtools} sort "
        "-@ {threads} "
        "-o {output.bam} -"


#Index BAM
rule samtools_index:
    input:
        bam = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam"

    output:
        bai = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam.bai"
    
    params:
        wd = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"],
        sample = "{reads}",
        samtools = config["DEPENDANCES"]["SAMTOOLS"]
    
    shell:
        "{params.samtools} index "
        "{input.bam}"



