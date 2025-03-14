rule bam_stat:
    input:
        bam = os.path.abspath(path_bam + name_genome + "/mapping/{reads}.sorted.bam"),
        csi = path_bam + name_genome + "/mapping/{reads}.sorted.bam.csi",

    output:        
        stats = path_qc + "BAM/" + name_genome +"/bam_stat/{reads}.stats",

    params:
        bam_stat = bam_stat

    shell:
        "{params.bam_stat} -i {input.bam} > {output.stats}"


rule read_gc:
    input:
        bam = os.path.abspath(path_bam + name_genome + "/mapping/{reads}.sorted.bam"),
        csi = path_bam + name_genome + "/mapping/{reads}.sorted.bam.csi",

    output:
        plot_pdf = path_qc + "BAM/" + name_genome + "/read_GC/{reads}.GC_plot.pdf",
        plot_r = path_qc + "BAM/" + name_genome + "/read_GC/{reads}.GC_plot.r",
        xls = path_qc + "BAM/" + name_genome + "/read_GC/{reads}.GC.xls"

    params:
        path_gen = path_qc + "BAM/" + name_genome + "/read_GC/{reads}",
        path = path_qc + "BAM/" + name_genome + "/read_GC/{reads}",
        read_GC  = read_GC

    shell:
        "mkdir -p {params.path_gen} && "
        "{params.read_GC} -i {input.bam} -o {params.path}"



rule samtools_stats:
    input:
        bam = os.path.abspath(path_bam + name_genome + "/mapping/{reads}.sorted.bam"),
        bai = path_bam + name_genome + "/mapping/{reads}.sorted.bam.bai",

    
    output:
        stats = path_qc + "BAM/" + name_genome + "/samtools/{reads}.stats"

    params:
        samtools = samtools

    threads:
        config["QUALITY_CONTROLE"]["THREADS"]

    shell:
        "{params.samtools} stats "
        "-@ {threads} {input.bam} > {output.stats}"



rule samtools_flagstat:
    input:
        bam = os.path.abspath(path_bam + name_genome + "/mapping/{reads}.sorted.bam"),
        bai = path_bam + name_genome + "/mapping/{reads}.sorted.bam.bai",
    
    output:
        flagstat = path_qc + "BAM/" + name_genome + "/samtools/{reads}.flagstat"

    params:
        samtools = samtools

    threads:
        config["QUALITY_CONTROLE"]["THREADS"]

    shell:
        "{params.samtools} flagstat "
        "-@ {threads} {input.bam} > {output.flagstat}"


# rule samtools_depth:
#     input:
#         bam = os.path.abspath(path_bam  + name_genome + "/{reads}.sorted.bam"),
#         bai = path_bam + name_genome + "/{reads}.sorted.bam.bai",
#         bed = path_bam + name_genome + "/{reads}_juncs.bed"

#     output:
#         depth = path_qc + "BAM/" + name_genome + "/samtools/{reads}.depth"

#     params:
#         samtools = samtools

#     shell:
#         "{params.samtools} depth "
#         "-b {input.bed} "
#         "{input.bam} "
#         "-o {output.depth}"

    
rule multiqc_bam:
    input:
        bam_stats = expand(path_qc + "BAM/" + name_genome +"/bam_stat/{reads}.stats",reads = all_samples),
        plot_pdf = expand(path_qc + "BAM/" + name_genome + "/read_GC/{reads}.GC_plot.pdf",reads = all_samples),
        plot_r = expand(path_qc + "BAM/" + name_genome + "/read_GC/{reads}.GC_plot.r",reads = all_samples),
        xls= expand(path_qc + "BAM/" + name_genome + "/read_GC/{reads}.GC.xls",reads = all_samples),
        samtools_stats = expand(path_qc + "BAM/" + name_genome + "/samtools/{reads}.stats",reads = all_samples),
        samtools_flagstats = expand(path_qc + "BAM/" + name_genome + "/samtools/{reads}.flagstat",reads = all_samples),
       #samtools_depth = expand(path_qc + "BAM/" + name_genome + "/samtools/{reads}.depth",reads = all_samples),
        log_final = expand(path_bam  + name_genome + "/mapping/log_star/{reads}_Log.final.out",reads = all_samples),
        log = expand(path_bam  + name_genome + "/mapping/log_star/{reads}_Log.out",reads = all_samples),
        log_progress = expand(path_bam  + name_genome + "/mapping/log_star/{reads}_Log.progress.out",reads = all_samples),
        tab = expand(path_bam  + name_genome + "/mapping/log_star/{reads}_SJ.out.tab",reads = all_samples),

    output:
        directory_data = directory(path_qc + "multiqc/BAM/" + name_genome + "/" + prefix +"_" + unique_id + "_data/"),
        html = path_qc + "multiqc/BAM/" + name_genome + "/" + prefix +"_" + unique_id + ".html"

    params:
        name = prefix + "_" + unique_id,
        multiqc = multiqc,
        path =  path_qc + "multiqc/BAM/" + name_genome + "/"

    threads:
        config["QUALITY_CONTROLE"]["THREADS"]
        
    shell:
        "{params.multiqc} "
        "{input.bam_stats} {input.plot_pdf} {input.plot_r} {input.xls} {input.samtools_stats} {input.samtools_flagstats} {input.log_final} {input.log} {input.log_progress} {input.tab} "
        "--filename {params.name} "
        "-o {params.path}"