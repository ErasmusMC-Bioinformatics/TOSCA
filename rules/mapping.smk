rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp(outputdir + "trimmed/{sample}_1_fastq.gz"),
        r2=temp(outputdir + "trimmed/{sample}_2_fastq.gz"),
        r1_unpaired=temp(outputdir + "trimmed/{sample}_1.unpaired.fastq.gz"),
        r2_unpaired=temp(outputdir + "trimmed/{sample}_2.unpaired.fastq.gz"),      
    params:
        **config["params"]["trimmomatic"]["pe"]
    threads:
        config["ncores"]    
    log:
        outputdir + "logs/trimmomatic/{sample}.log"
    benchmark:
        outputdir + "benchmarks/trimmomatic/{sample}.trim.benchmark.txt"
    wrapper:
        "v1.23.3/bio/trimmomatic/pe"

rule map_reads:
    input:
        fq1=get_trimmed_reads,
        # path to STAR reference genome index
        idx=rules.star_index.output
    output:
        # see STAR manual for additional output files
        log=outputdir + "logs/star/{sample}/Log.out",
        sj=outputdir + "star/pe/{sample}/SJ.out.tab",
	    aln = temp(outputdir + "mapped/{sample}.sorted.bam")
    log:
        outputdir + "logs/star/{sample}.log"
    benchmark:
        outputdir + "benchmarks/star/{sample}.star.benchmark.txt"
    params:
        # optional parameters
        twopass="--twopassMode Basic"
    threads: config["ncores"]
    wrapper:
        "v1.23.3/bio/star/align"

rule samtools_sort:
    input:
        outputdir + "mapped/{sample}.bam",
    output:
        outputdir + temp("mapped/{sample}.sorted.bam"),
    log:
        "{sample}.log",
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        "v1.23.0/bio/samtools/sort"        

rule add_read_groups:
    input:
        outputdir + "mapped/{sample}.sorted.bam",
    output:
        temp(outputdir + "mapped/fixed-rg/{sample}.bam"),
    log:
        outputdir + "logs/picard/replace_rg/{sample}.log",
    params:
        extra="--RGLB lib1 --RGPL ILLUMINA --RGPU {sample} --RGSM {sample}",
    resources:
        mem_mb=1024,
    wrapper:
        "master/bio/picard/addorreplacereadgroups"


rule mark_duplicates:
    input:
        bams=outputdir + "mapped/fixed-rg/{sample}.bam"
    output:
        bam=temp(outputdir + "dedup/{sample}.bam"),
        metrics=outputdir + "qc/dedup/{sample}.metrics.txt"
    log:
        outputdir + "logs/picard/dedup/{sample}.log"
    benchmark:
        outputdir + "benchmarks/picard/{sample}.picard.benchmark.txt"
    params:
        extra="--REMOVE_DUPLICATES=true"
    resources:
        mem_mb=4096    
    wrapper:
        "v1.23.3/bio/picard/markduplicates"

rule splitncigarreads:
    input:
        bam=outputdir + "dedup/{sample}.bam",
        ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
    output:
        temp(outputdir + "split/{sample}.bam"),
    log:
        outputdir + "logs/gatk/splitNCIGARreads/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "master/bio/gatk/splitncigarreads"        
    

rule gatk_baserecalibrator:
    input:
        bam=outputdir + "split/{sample}.bam",
        ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
        dict=expand("resources/reference_genome/{ref}/{species}.dict",ref=config["ref"]["build"],species=config["ref"]["species"]),
        known=expand("resources/database/{ref}/variation.noiupac.vcf.gz",ref=config["ref"]["build"]),
        tbi=expand("resources/database/{ref}/variation.noiupac.vcf.gz.tbi",ref=config["ref"]["build"])
    output:
        recal_table=outputdir + "recal/{sample}.table"
    log:
        outputdir + "logs/gatk/bqsr/{sample}.log"
    benchmark:
        outputdir + "benchmarks/gatk/bqsr/{sample}.gatk.bqsr.benchmark.txt"
    wrapper:
        "0.78.0/bio/gatk/baserecalibrator"    

rule gatk_applybqsr:
    input:
        bam=outputdir + "split/{sample}.bam",
        ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
        dict=expand("resources/reference_genome/{ref}/{species}.dict",ref=config["ref"]["build"],species=config["ref"]["species"]),
        recal_table=outputdir + "recal/{sample}.table"
    output:
        bam=outputdir + "recal/{sample}.bam"
    log:
        outputdir + "logs/gatk/gatk_applybqsr/{sample}.log"
    benchmark:
        outputdir + "benchmarks/gatk/gatk_applybqsr/{sample}.gatk.applybqsr.benchmark.txt"
    wrapper:
        "0.78.0/bio/gatk/applybqsr"
