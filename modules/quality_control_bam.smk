rule bam_stat:
    input:
        bam = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam"

    output:        
        stats = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/bam_stat/{reads}.stats"

    shell:
        "bam_stat.py -i {input.bam} > {output.stats}"


rule read_gc:
    input:
        bam = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam"
    
    output:
        plot_pdf = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/read_GC/{reads}.GC_plot.pdf",
        plot_r = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/read_GC/{reads}.GC_plot.r",
        xls= config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/read_GC/{reads}.GC.xls"

    params:
        path = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/read_GC/{reads}"

    shell:
        "read_GC.py -i {input.bam} -o {params.path}"

# rule kmerexplor:
#     input:
#         read = expand(config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/reads1/{reads}_1.fastq.gz"),
       

#     output:

#     params:
#         fastqc = config["DEPENDANCES"]["KMEREXPLOR"]



#     threads:
#         config["PARAMS"]["KMEREXPLOR"]["THREADS"]

#     shell:


# rule qualimap:
#     input:

#     output:

#     params:

#     threads:

#     shell:

rule samtools_stats:
    input:
        bam = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam"
    
    output:
        stats = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/samtools/{reads}.stats"

    params:
        samtools = config["DEPENDANCES"]["SAMTOOLS"]

    threads:
        config["PARAMS"]["SAMTOOLS"]["THREADS"]

    shell:
        "{params.samtools} stats "
        "-@ {threads} {input.bam} > {output.stats}"


rule samtools_flagstat:
    input:
        bam = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam"
    
    output:
        flagstat = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/samtools/{reads}.flagstat"

    params:
        samtools = config["DEPENDANCES"]["SAMTOOLS"]

    threads:
        config["PARAMS"]["SAMTOOLS"]["THREADS"]

    shell:
        "{params.samtools} flagstat "
        "-@ {threads} {input.bam} > {output.flagstat}"



rule multiqc:
    input:
        html1 = expand(config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc/{reads}_1_fastqc.html",reads = all_samples),
        zip1 = expand(config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc/{reads}_1_fastqc.zip", reads = all_samples),
        html2 = expand(config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc/{reads}_2_fastqc.html",reads = all_samples),
        zip2 = expand(config["DATA_INPUT"]["WORKING_DIRECTORY"] + "/2-processed_data/quality_control/fastqc/{reads}_2_fastqc.zip", reads = all_samples),
        bam_stats = expand(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/bam_stat/{reads}.stats",reads = all_samples),
        plot_pdf = expand(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/read_GC/{reads}.GC_plot.pdf",reads = all_samples),
        plot_r = expand(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/read_GC/{reads}.GC_plot.r",reads = all_samples),
        xls= expand(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/read_GC/{reads}.GC.xls",reads = all_samples),
        samtools_stats = expand(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/samtools/{reads}.stats",reads = all_samples),
        flagstats = expand(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/samtools/{reads}.flagstat",reads = all_samples)

    output:
        directory_data = directory(config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/multiqc/" + config["PARAMS"]["GENERAL"]["PREFIX"] +"_data/"),
        html = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/multiqc/" + config["PARAMS"]["GENERAL"]["PREFIX"]  + ".html"
    
    params:
        name = config["PARAMS"]["GENERAL"]["PREFIX"],
        multiqc = config["DEPENDANCES"]["MULTIQC"],
        path_raw = config["DATA_INPUT"]["WORKING_DIRECTORY"],
        path_analysis = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"],
        path = config["PARAMS"]["GENERAL"]["WORKING_DIRECTORY"] + config["PARAMS"]["GENERAL"]["PREFIX"] + "/quality_control/multiqc/"

    threads:
        config["PARAMS"]["MULTIQC"]["THREADS"]
        
    shell:
        "{params.multiqc} "
        "{params.path_raw} {params.path_analysis} "
        "--filename {params.name} "
        "-o {params.path}"


