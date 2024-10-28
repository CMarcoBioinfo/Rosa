rule bam_stat:
    input:
        bam = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam"

    output:        
        stats = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/bam_stat/{reads}.stats"

    shell:
        "bam_stat.py -i {input.bam} > {output.stats}"


rule read_gc:
    input:
        bam = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam"
    
    output:
        plot_pdf = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/read_GC/{reads}.GC_plot.pdf",
        plot_r = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/read_GC/{reads}.GC_plot.r",
        xls= config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/read_GC/{reads}.GC.xls"

    params:
        path = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/read_GC/{reads}"

    shell:
        "read_GC.py -i {input.bam} -o {params.path}"

# rule kmerexplor:
#     input:
#         read = expand(config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/samples/reads1/{reads}_1.fastq.gz"),
       

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
        bam = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam",
        bai = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam.bai"

    
    output:
        stats = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/samtools/{reads}.stats"

    params:
        samtools = config["DEPENDANCES"]["MAPPING"]["SAMTOOLS"]

    threads:
        config["PARAMS"]["SAMTOOLS"]["THREADS"]

    shell:
        "{params.samtools} stats "
        "-@ {threads} {input.bam} > {output.stats}"



rule samtools_flagstat:
    input:
        bam = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam",
        bai = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam.bai"
    
    output:
        flagstat = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/samtools/{reads}.flagstat"

    params:
        samtools = config["DEPENDANCES"]["MAPPING"]["SAMTOOLS"]

    threads:
        config["PARAMS"]["SAMTOOLS"]["THREADS"]

    shell:
        "{params.samtools} flagstat "
        "-@ {threads} {input.bam} > {output.flagstat}"


rule samtools_depth:
    input:
        bam = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam",
        bai = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam.bai",
        bed = (config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/annotation/" + config["GENERAL"]["DATA_INPUTS"]["ANNOTATION"]).rsplit(".",1)[0] + ".bed"


    output:
        depth = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/samtools/{reads}.depth"

    params:
        samtools = config["DEPENDANCES"]["MAPPING"]["SAMTOOLS"]

    shell:
        "{params.samtools} depth "
        "-b {input.bed} "
        "{input.bam} "
        "-o {output.depth}"


# rule samtools_coverage:
#     input:
#         bam = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam",
#         bai = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/1-mapping/{reads}.sorted.bam.bai",
#         bed3 = (config["GENERAL"]["DATA_INPUTS"]["WORKING_DIRECTORY"] + "/1-raw_data/annotation/" + config["GENERAL"]["DATA_INPUTS"]["ANNOTATION"]).rsplit(".",1)[0] + ".bed3"


#     output:
#         cov = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/samtools/{reads}.cov"

#     params:
#         samtools = config["DEPENDANCES"]["MAPPING"]["SAMTOOLS"]

#     shell:
#         "{params.samtools} coverage "
#         "-b {input.bed3} "
#         "{input.bam} "
#         "-o {output.cov}"


# rule qualimap:
#     input:

#     output:

#     params:


#     threads:

#     shell:


    
rule multiqc_bam:
    input:
        bam_stats = expand(config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/bam_stat/{reads}.stats",reads = all_samples),
        plot_pdf = expand(config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/read_GC/{reads}.GC_plot.pdf",reads = all_samples),
        plot_r = expand(config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/read_GC/{reads}.GC_plot.r",reads = all_samples),
        xls= expand(config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/read_GC/{reads}.GC.xls",reads = all_samples),
        samtools_stats = expand(config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/samtools/{reads}.stats",reads = all_samples),
        samtools_flagstats = expand(config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/samtools/{reads}.flagstat",reads = all_samples),
        samtools_depth = expand(config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/samtools/{reads}.depth",reads = all_samples)

    output:
        directory_data = directory(config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/multiqc/" + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] +"_bam_" + unique_id + "_data/"),
        html = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/multiqc/" + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"]  + "_bam_" + unique_id + ".html"
    
    params:
        name = config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "_bam_" + unique_id,
        multiqc = config["DEPENDANCES"]["QC"]["MULTIQC"],
        path_analysis = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"],
        path = config["GENERAL"]["DATA_OUTPUTS"]["WORKING_DIRECTORY"] + config["GENERAL"]["DATA_OUTPUTS"]["PREFIX"] + "/quality_control/multiqc/"

    threads:
        config["PARAMS"]["MULTIQC"]["THREADS"]
        
    shell:
        "{params.multiqc} "
        "{input.bam_stats} {input.plot_pdf} {input.plot_r} {input.xls} {input.samtools_stats} {input.samtools_flagstats} {input.samtools_depth} "
        "--filename {params.name} "
        "-o {params.path}"