# index reference faidx

# samtools faidx 1-raw_data/reference/Homo_sapiens.GRCh37.dna.primary_assembly.chr.fa

# index VCF: tabix



# Gatk IndexFeatureFile

# ~/Bureau/gatk-4.6.1.0/gatk IndexFeatureFile -I 1-raw_data/annotation/00-All.chr.vcf.gz

# Ajout groupe de lecture (peut etre fait pendant le mapping) + changer le fichier de configuration

# gatk CreateDictionary 

# ~/Bureau/gatk-4.6.1.0/gatk CreateSequenceDictionary -R /home/m_corentin/Rosa/data/1-raw_data/reference/Homo_sapiens.GRCh37.dna.primary_assembly.chr.fa

# gatk SpliceNcigarReads

# ~/Bureau/gatk-4.6.1.0/gatk SplitNCigarReads -R 1-raw_data/reference/Homo_sapiens.GRCh37.dna.primary_assembly.chr.fa -I 2-processed_data/trimmed/BAM/spliceLauncher_STAR/Homo_sapiens.GRCh37.dna.primary_assembly.chr/2006161710.100bp.sorted.bam -O test/2006161710.100bp.sorted.split.bam


# BaseRecalibrator

# ~/Bureau/gatk-4.6.1.0/gatk BaseRecalibrator -I test/2006161710.100bp.sorted.split.rg.bam -R 1-raw_data/reference/Homo_sapiens.GRCh37.dna.primary_assembly.chr.fa --bqsr-recal-file test/2006161710.100bp.sorted.split.recal_data.table -O test/2006161710.100bp.sorted.split.rg.recalibrated.bam

# ApplyBQRS

# ~/Bureau/gatk-4.6.1.0/gatk ApplyBQSR -I test/2006161710.100bp.sorted.split.rg.bam -R 1-raw_data/reference/Homo_sapiens.GRCh37.dna.primary_assembly.chr.fa --bqsr-recal-file test/2006161710.100bp.sorted.split.recal_data.table -O test/2006161710.100bp.sorted.split.rg.recalibrated.bam

ruleorder: star_align > gatk_SplitNCigarReads

if use_trimming:
    sequencing_date = lambda wildcards: dict_csv.loc[dict_csv['id'] == str(wildcards.reads).rsplit('.',1)[0], 'sequencing_date'].values[0]
    technologie = lambda wildcards: dict_csv.loc[dict_csv['id'] == str(wildcards.reads).rsplit('.',1)[0], 'technologie'].values[0]
else : 
    sequencing_date = lambda wildcards: dict_csv.loc[dict_csv['id'] == str(wildcards.reads), 'sequencing_date'].values[0]
    technologie = lambda wildcards: dict_csv.loc[dict_csv['id'] == str(wildcards.reads), 'technologie'].values[0]



rule gatk_AddOrReplaceReadGroups:
    input:
        bam = os.path.abspath(path_bam + name_genome + "/mapping/{reads}.sorted.bam"),
        csi = path_bam + name_genome + "/mapping/{reads}.sorted.bam.csi",

    output:
        bam_ReadsGroups = temp(os.path.abspath(path_bam  + name_genome + "/GATK/{reads}.ReadsGroups.sorted.bam")),
        bai = temp(path_bam + name_genome + "/GATK/{reads}.ReadsGroups.sorted.bam.bai"),
        csi = temp(path_bam + name_genome + "/GATK/{reads}.ReadsGroups.sorted.bam.csi")

    params:
        bam_ReadsGroups = os.path.abspath(path_bam  + name_genome + "/GATK/{reads}.ReadsGroups.bam"),
        gatk = gatk,
        samtools = samtools,
        sequencing_date = sequencing_date,
        technologie = technologie

    shell: 
        "{params.gatk} AddOrReplaceReadGroups "
        "--INPUT {input.bam} "
        "--OUTPUT {params.bam_ReadsGroups} "
        "--RGID {params.sequencing_date} "
        "--RGLB lib_{params.sequencing_date} "
        "--RGPL {params.technologie} "
        "--RGPU 1234 "
        "--RGSM {wildcards.reads} && "
        "{params.samtools} sort "
        "-@ {threads} "
        "-o {output.bam_ReadsGroups} "
        "--output-fmt BAM --write-index "
        "{params.bam_ReadsGroups} && "
        "{params.samtools} index "
        "-@ {threads} "
        "-b {output.bam_ReadsGroups} && "
        "rm {params.bam_ReadsGroups}"

rule bcftools_index_vcf_known_sites:
    input:
        vcf_known_sites = vcf_known_sites

    output:
        tbi = vcf_known_sites + ".tbi"

    params:
        bcftools = bcftools

    shell:
        "{params.bcftools} index {input.vcf_known_sites}"


rule samtools_faidx:
    input: 
        reference = genome

    output:
        fai = genome + ".fai"
    
    params:
        samtools = samtools

    threads:
        config["MAPPING"]["INDEX"]["THREADS"]

    shell:
        "{params.samtools} faidx "
        "-@ {threads}"
        "{input.reference}"


rule gatk_CreateSequenceDictionary:
    input:
        reference = genome,
        fai = genome + ".fai"

    output:
        dict_file = genome.rsplit(".",1)[0] + ".dict" 

    params:
        gatk = gatk

    threads:
        config["MAPPING"]["INDEX"]["THREADS"]
    
    shell:
        """if [ -f {output.dict_file} ]; then
            rm {output.dict_file} 
        fi
        {params.gatk}  --java-options "-XX:ParallelGCThreads={threads}" CreateSequenceDictionary -R {input.reference}"""


# rule gatk_SplitNCigarReads:
#     input:
#         reference = genome,
#         dict_file = genome.rsplit(".",1)[0] + ".dict" ,
#         bam = os.path.abspath(path_bam  + name_genome + "/{reads}.ReadsGroups.bam")
            
#     output:
#         bam_split =  os.path.abspath(path_bam  + name_genome + "/{reads}.split.bam")

#     params:
#         gatk = gatk

#     threads:
#         config["MAPPING"]["RECALIBRATE"]["THREADS"]

#     shell:
#         "{params.gatk} "
#         "SplitNCigarReads "
#         "-R {input.reference} "
#         "-I {input.bam} "
#         "-O {output.bam_split}"


# rule gatk_MarkDuplicates :
#     input:
#         bam = os.path.abspath(path_bam  + name_genome + "/GATK/{reads}.ReadsGroups.sorted.bam"),
#         bai = path_bam + name_genome + "/GATK/{reads}.ReadsGroups.sorted.bam.bai",
#         csi = path_bam + name_genome + "/GATK/{reads}.ReadsGroups.sorted.bam.csi"

#     output:
#         bam_marked_sorted = temp(os.path.abspath(path_bam + name_genome + "/GATK/{reads}.marked_duplicates.sorted.bam")),
#         bai = temp(path_bam + name_genome + "/GATK/{reads}.marked_duplicates.sorted.bam.bai"),
#         csi = temp(path_bam + name_genome + "/GATK/{reads}.marked_duplicates.sorted.bam.csi"),
#         metrics = path_bam + name_genome + "/GATK/{reads}.marked_dup_metrics.txt"

#     params:
#         gatk = gatk,
#         samtools = samtools,
#         tmp_dir = path_bam + name_genome + "/GATK/{reads}_MarkDuplicates"

#     threads:
#         config["GATK"]["RECALIBRATE"]["THREADS"]

#     shell:
#         "mkdir -p {params.tmp_dir} && "
#         "samtools idxstats {input.bam} | cut -f1 | grep -v '*' | parallel -j {threads} '"
#         "mkdir -p {params.tmp_dir}/tmp_dir_{{}} && "
#         "{params.samtools} view -b {input.bam} {{}} > {params.tmp_dir}/{{}}.bam && "
#         "{params.samtools} index {params.tmp_dir}/{{}}.bam && "
#         "{params.gatk} MarkDuplicates "
#         "-I {params.tmp_dir}/{{}}.bam "
#         "-O {params.tmp_dir}/marked_reads_{{}}.bam "
#         "-M {params.tmp_dir}/marked_reads_{{}}.metrics.txt "
#         "--TMP_DIR {params.tmp_dir}/tmp_dir_{{}}' && "
#         "{params.samtools} merge "
#         "-@ {threads} -f {params.tmp_dir}/merge.bam "
#         "{params.tmp_dir}/marked_reads_*.bam && "
#         "{params.samtools} sort "
#         "-@ {threads} "
#         "-o {output.bam_marked_sorted} "
#         "--output-fmt BAM --write-index "
#         "{params.tmp_dir}/merge.bam && "
#         "{params.samtools} index "
#         "-@ {threads} "
#         "-b {output.bam_marked_sorted} && "
#         "cat {params.tmp_dir}/marked_reads_*.metrics.txt > {output.metrics} && "
#         "rm -r {params.tmp_dir}"


rule gatk_MarkDuplicates :
    input:
        bam = os.path.abspath(path_bam  + name_genome + "/GATK/{reads}.ReadsGroups.sorted.bam"),
        bai = path_bam + name_genome + "/GATK/{reads}.ReadsGroups.sorted.bam.bai",
        csi = path_bam + name_genome + "/GATK/{reads}.ReadsGroups.sorted.bam.csi"

    output:
        bam_marked_sorted = temp(os.path.abspath(path_bam + name_genome + "/GATK/{reads}.marked_duplicates.sorted.bam")),
        bai = temp(path_bam + name_genome + "/GATK/{reads}.marked_duplicates.sorted.bam.bai"),
        csi = temp(path_bam + name_genome + "/GATK/{reads}.marked_duplicates.sorted.bam.csi"),
        metrics = path_bam + name_genome + "/GATK/{reads}.marked_dup_metrics.txt"

    params:
        gatk = gatk,
        samtools = samtools,
        tmp_dir = path_bam + name_genome + "/GATK/{reads}_MarkDuplicates",
        bam_marked =  path_bam + name_genome + "/GATK/{reads}_MarkDuplicates/{reads}.marked_duplicates.bam",


    shell:
        "mkdir -p {params.tmp_dir} && "
        "{params.gatk} MarkDuplicates "
        "-I {input.bam} "
        "-O {params.bam_marked} "
        "-M {output.metrics} "
        "--TMP_DIR {params.tmp_dir} && "
        "{params.samtools} sort "
        "-o {output.bam_marked_sorted} "
        "--output-fmt BAM --write-index "
        "{params.bam_marked} && "
        "{params.samtools} index "
        "-b {output.bam_marked_sorted} && "
        "rm -r {params.tmp_dir}"


rule gatk_SplitNCigarReads:
    input:
        reference = genome,
        dict_file = genome.rsplit(".",1)[0] + ".dict" ,
        bam = os.path.abspath(path_bam + name_genome + "/GATK/{reads}.marked_duplicates.sorted.bam"),
        bai = path_bam + name_genome + "/GATK/{reads}.marked_duplicates.sorted.bam.bai",
        csi = path_bam + name_genome + "/GATK/{reads}.marked_duplicates.sorted.bam.csi",
            
    output:
        bam_split_sorted =  temp(os.path.abspath(path_bam  + name_genome + "/GATK/{reads}.split.sorted.bam")),
        bai =  temp(path_bam  + name_genome + "/GATK/{reads}.split.sorted.bam.bai"),
        csi =  temp(path_bam  + name_genome + "/GATK/{reads}.split.sorted.bam.csi")


    params:
        gatk = gatk,
        samtools = samtools,
        tmp_dir = path_bam  + name_genome + "/GATK/{reads}_SplitNCigarReads"

    threads:
        config["GATK"]["RECALIBRATE"]["THREADS"]

    shell:
        "mkdir -p {params.tmp_dir} && "
        "samtools idxstats {input.bam} | cut -f1 | grep -v '*' | parallel -j {threads} '"
        "mkdir -p {params.tmp_dir}/tmp_dir_{{}} && "
        "{params.gatk} SplitNCigarReads "
        "-R {input.reference} "
        "-I {input.bam} "
        "-O {params.tmp_dir}/split_reads_{{}}.bam "
        "-L {{}} "
        "--tmp-dir {params.tmp_dir}/tmp_dir_{{}}' && "
        "{params.samtools} merge "
        "-@ {threads} -f {params.tmp_dir}/merge.bam "
        "{params.tmp_dir}/split_reads_*.bam && "
        "{params.samtools} sort "
        "-@ {threads} "
        "-o {output.bam_split_sorted} "
        "--output-fmt BAM --write-index "
        "{params.tmp_dir}/merge.bam && "
        "{params.samtools} index "
        "-@ {threads} "
        "-b {output.bam_split_sorted} && "
        "rm -r {params.tmp_dir}"



rule gatk_BaseRecalibrator:
    input:
        reference = genome,
        bam =  os.path.abspath(path_bam  + name_genome + "/GATK/{reads}.split.sorted.bam"),
        bai =  path_bam  + name_genome + "/GATK/{reads}.split.sorted.bam.bai",
        csi =  path_bam  + name_genome + "/GATK/{reads}.split.sorted.bam.csi",
        vcf_known_sites = vcf_known_sites

    output:
        table = temp(os.path.abspath(path_bam  + name_genome + "/GATK/{reads}.recal.table")),

    params:
        gatk = gatk,

    threads:
        config["GATK"]["RECALIBRATE"]["THREADS"]


    shell:
        "{params.gatk} BaseRecalibrator "
        "-R {input.reference} "
        "-I {input.bam} "
        "--known-sites {input.vcf_known_sites} "
        "-O {output.table}"


rule gatk_ApplyBQSR:
    input:
        reference = genome,
        bam =  os.path.abspath(path_bam  + name_genome + "/GATK/{reads}.split.sorted.bam"),
        bai =  path_bam  + name_genome + "/GATK/{reads}.split.sorted.bam.bai",
        table = os.path.abspath(path_bam  + name_genome + "/GATK/{reads}.recal.table")

    output:
        bam_recalibrated = os.path.abspath(path_bam  + name_genome + "/GATK/{reads}.recalibrated.bam"),
        bai = path_bam  + name_genome + "/GATK/{reads}.recalibrated.bam.bai",
        csi = path_bam  + name_genome + "/GATK/{reads}.recalibrated.bam.csi"


    params:
        gatk = gatk,
        samtools = samtools,
        tmp_dir = path_bam  + name_genome + "/GATK/{reads}_ApllyBQSR"

    threads:
        config["GATK"]["RECALIBRATE"]["THREADS"]

    shell:
        "mkdir -p {params.tmp_dir} && "
        "samtools idxstats {input.bam} | cut -f1 | grep -v '*' | parallel -j {threads} '"
        "mkdir -p {params.tmp_dir}/tmp_dir_{{}} && "
        "{params.gatk} ApplyBQSR "
        "-R {input.reference} "
        "-I {input.bam} "
        "--bqsr-recal-file {input.table} "
        "-O {params.tmp_dir}/recalibrated_{{}}.bam "
        "-L {{}} "
        "--tmp-dir {params.tmp_dir}/tmp_dir_{{}}' && "
         "{params.samtools} merge "
        "-@ {threads} -f {params.tmp_dir}/merge.recalibrated.bam "
        "{params.tmp_dir}/recalibrated_*.bam && "
        "{params.samtools} sort "
        "-@ {threads} "
        "-o {output.bam_recalibrated} "
        "--output-fmt BAM --write-index "
        "{params.tmp_dir}/merge.recalibrated.bam && "
        "{params.samtools} index "
        "-@ {threads} "
        "-b {output.bam_recalibrated} && "
        "rm -r {params.tmp_dir}"



# rule samtools_sort_recalibrated_bam:
#     input:
#         bam_recalibrated = os.path.abspath(path_bam  + name_genome + "/{reads}.recalibrated.bam")

#     output:
#         sort_bam = os.path.abspath(path_bam  + name_genome + "/{reads}.recalibrated.sorted.bam"),
#         csi = path_bam  + name_genome + "/{reads}.recalibrated.sorted.bam.csi",

#     params:
#         samtools = samtools

#     threads:
#         config["MAPPING"]["ALIGN"]["THREADS"]

#     shell:
#         "{params.samtools} sort "
#         "-o {output.sort_bam} "
#         "--output-fmt BAM "
#         "--write-index "
#         "--threads {threads} "
#         "{input.bam_recalibrated}"

# rule samtools_index:
#     input:
#         sort_bam = os.path.abspath(path_bam  + name_genome + "/{reads}.recalibrated.sorted.bam"),

#     output:
#         bai = path_bam  + name_genome + "/{reads}.recalibrated.sorted.bam.bai",

#     params:
#         samtools = samtools

#     threads:
#         config["MAPPING"]["ALIGN"]["THREADS"]

#     shell:
#         "{params.samtools} index "
#         "-@ {threads}"
#         "-b {input.sort_bam}"