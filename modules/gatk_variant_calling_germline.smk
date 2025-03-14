def parse_dict_file(dict_file):
    seq_lengths = {}

    with open(dict_file, "r") as file:
        for line in file:
            items = line.strip().split()
            SN_item = [item for item in items if item.startswith("SN:")]
            LN_item = [item for item in items if item.startswith("LN:")]

            if SN_item and LN_item:
                chr = SN_item[0].split("SN:")[1]
                length = int(LN_item[0].split("LN:")[1])
                seq_lengths[chr] = length
    
    return seq_lengths

def max_chromosome_length(seq_lengths):
    return max(seq_lengths.values())

def threads_parallel (threads):
    threads_use = int(threads / 4)
    if threads_use < 1:
        threads_use = 1
    return threads_use

dict_file = genome.rsplit(".",1)[0] + ".dict" 
seq_lengths = parse_dict_file(dict_file)
max_chromosome_length = max_chromosome_length(seq_lengths)
memory_needed_mb_HaplotypeCaller = estimate_memory(max_chromosome_length, 11)


rule gatk_HaplotypeCaller:
    input:
        reference = genome,
        bam = os.path.abspath(path_bam  + name_genome + "/GATK/{reads}.recalibrated.bam"),
        bai = path_bam  + name_genome + "/GATK/{reads}.recalibrated.bam.bai",
        csi = path_bam  + name_genome + "/GATK/{reads}.recalibrated.bam.csi"

    output:
        vcf = path_results + "/GATK/variant_calling/{reads}.HaplotypeCaller.vcf.gz",
        tbi = path_results + "/GATK/variant_calling/{reads}.HaplotypeCaller.vcf.gz.tbi",

    params:
        gatk = gatk,
        tmp_dir = path_results + "/GATK/variant_calling/{reads}_HaplotypeCaller",


    threads:
        lambda wildcards: calculate_optimal_threads(mem_total_mb, memory_needed_mb_HaplotypeCaller, max_cores, int(config["GATK"]["VARIANT_CALLING"]["THREADS"]))

    # resources:
    #     mem_mb = lambda wildcards,threads: calculate_mem_per_rule(mem_total_mb, memory_needed_mb_HaplotypeCaller, threads)  

    run:
        threads_use = threads_parallel(threads)
        shell("mkdir -p {params.tmp_dir} && "
        "samtools idxstats {input.bam} | cut -f1 | grep -v '*' | parallel -j {threads_use} '"
        "mkdir -p {params.tmp_dir}/tmp_dir_{{}} && "
        "{params.gatk} HaplotypeCaller "
        "-R {input.reference} "
        "-I {input.bam} "
        "-O {params.tmp_dir}/HaplotypeCaller_{{}}.vcf.gz "
        "-L {{}} "
        "--tmp-dir {params.tmp_dir}/tmp_dir_{{}}'")
        shell("{params.gatk} MergeVcfs "
        "$(for vcf in {params.tmp_dir}/HaplotypeCaller_*.vcf.gz; do echo -I $vcf; done) "
        "-O {output.vcf}")
        shell("rm -r {params.tmp_dir}")


LowDP = config["GATK"]["VARIANT_CALLING"].get("LOW_DP")
if not LowDP:
    LowDP = float(20.0)

LowMQ = config["GATK"]["VARIANT_CALLING"].get("LOW_MQ")
if not LowMQ:
    LowMQ = float(55.0)

HighFS = config["GATK"]["VARIANT_CALLING"].get("HIGH_FS")
if not HighFS:
    HighFS = float(30.0)

LowQD = config["GATK"]["VARIANT_CALLING"].get("LOW_QD")
if not LowQD:
    LowQD = float(5.0)

HighSOR = config["GATK"]["VARIANT_CALLING"].get("HIGH_SOR")
if not HighSOR:
    HighSOR = float(3.0)

rule gatk_VariantFiltration:
    input:
        reference = genome,
        vcf = path_results + "/GATK/variant_calling/{reads}.HaplotypeCaller.vcf.gz",
        tbi = path_results + "/GATK/variant_calling/{reads}.HaplotypeCaller.vcf.gz.tbi",

    output:
        vcf = path_results + "/GATK/variant_calling/{reads}.VariantFiltration.vcf.gz",
        tbi = path_results + "/GATK/variant_calling/{reads}.VariantFiltration.vcf.gz.tbi"


    params:
        gatk = gatk,
        LowDP = LowDP,
        LowMQ = LowMQ,
        HighFS = HighFS,
        LowQD = LowQD,
        HighSOR = HighSOR

    shell:
        "{params.gatk} VariantFiltration "
        "-R {input.reference} "
        "-V {input.vcf} "
        "-O {output.vcf} "
        '--filter-expression "DP < {params.LowDP}" --filter-name "LowDP" '
        '--filter-expression "MQ < {params.LowMQ}" --filter-name "LowMQ" '
        '--filter-expression "FS > {params.HighFS}" --filter-name "HighFS" '
        '--filter-expression "QD < {params.LowQD}" --filter-name "LowQD" '
        '--filter-expression "SOR > {params.HighSOR}" --filter-name "HighSOR"'



rule pass_filter_variants:
    input:        
        vcf = path_results + "/GATK/variant_calling/{reads}.VariantFiltration.vcf.gz",
        tbi = path_results + "/GATK/variant_calling/{reads}.VariantFiltration.vcf.gz.tbi",

    output:
        vcf = path_results + "/GATK/variant_calling/{reads}.filtered.vcf.gz",
        tbi = path_results + "/GATK/variant_calling/{reads}.filtered.vcf.gz.tbi",

    params:
        bcftools = bcftools,
        bgzip = bgzip
    
    shell:
        "{params.bcftools} view "
        '-f "PASS" {input.vcf} | '
        "{params.bgzip} -c > {output.vcf} && "
        "{params.bcftools} index -t {output.vcf}"
