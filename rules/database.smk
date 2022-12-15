##### Download database #####

rule download_database:
    output:
        Gen1K = expand("resources/database/{ref}/1000G-phase_3.vcf.gz",ref=config["ref"]["build"]),
        ClinVar = expand("resources/database/{ref}/ClinVar.vcf.gz",ref=config["ref"]["build"]),
	ClinVar_tbi = expand("resources/database/{ref}/ClinVar.vcf.gz.tbi",ref=config["ref"]["build"]),
        COSMIC = expand("resources/database/{ref}/COSMIC.vcf.gz",ref=config["ref"]["build"])
    params:
        build=config["ref"]["build"],
        Cosmic = config["database_url"]["GRCh38"]["somatic"]["Cosmic"] if config["ref"]["build"]=='GRCh38' else config["database_url"]["GRCh37"]["somatic"]["Cosmic"],
	Gen1K = config["database_url"]["GRCh38"]["germline"]["1000G"] if config["ref"]["build"]=='GRCh38' else config["database_url"]["GRCh37"]["germline"]["1000G"],
	Clinvar = config["database_url"]["GRCh38"]["germline"]["ClinVar"] if config["ref"]["build"]=='GRCh38' else config["database_url"]["GRCh37"]["germline"]["ClinVar"],
    shell:
        "curl -k -L  '{params.Cosmic}' > {output.COSMIC}; "
        "curl -k -L  '{params.Gen1K}' > {output.Gen1K}; "
        "curl -k -L  '{params.Clinvar}' > {output.ClinVar}; "
        "curl -k -L  '{params.Clinvar}.tbi' > {output.ClinVar_tbi}; "

rule download_mapping:
    output:
        expand("resources/database/{ref}/raw_mappability_{len}.bw",ref=config["ref"]["build"],len=config["kmer"])
    params:
        build=config["ref"]["build"],
      	mappability = config["mappability"]["GRCh38"] if config["ref"]["build"]=='GRCh38' else config["mappability"]["GRCh37"],
        kmer = config["kmer"]
    cache: True
    run:
        if config["ref"]["build"]=='GRCh38': 
          shell("curl -L  '{params.mappability}_{params.kmer}.bw' > resources/database/{params.build}/raw_mappability_{params.kmer}.bw; ")
        else:
          shell("curl -L  '{params.mappability}{params.kmer}mer.bigWig' > resources/database/{params.build}/raw_mappability_{params.kmer}.bw; ")
        	
rule ESP6500SI:
    input:
         ESP = config["database_url"]["GRCh38"]["germline"]["ESP"],
         script = "scripts/ESP6500SI.R"
    output:
        expand("resources/database/{ref}/ESP6500SI.vcf.gz",ref=config["ref"]["build"])
    params:
        build=config["ref"]["build"]      
    cache: True
    conda:
        "../envs/r4.yaml"
    shell:
        "tar -xzvf {input.ESP} -C resources/database/{params.build}; "
        "Rscript --vanilla {input.script} {params.build}; "
        "vcf-concat resources/database/{params.build}/*.vcf | sort -k1,1V -k2,2n | bgzip -c > resources/database/{params.build}/ESP6500SI.vcf.gz; "
        "rm resources/database/{params.build}/*.vcf; "
	
rule tabix_Cosmic:
    input:
        expand("resources/database/{ref}/COSMIC.vcf.gz",ref=config["ref"]["build"])
    output:
        expand("resources/database/{ref}/COSMIC.vcf.gz.tbi",ref=config["ref"]["build"])
    params:
        "-p vcf"
    cache: True
    wrapper:
        "master/bio/tabix/index" #no index and 0.78.0 

rule tabix_1000G:
    input:
        expand("resources/database/{ref}/1000G-phase_3.vcf.gz",ref=config["ref"]["build"])
    output:
       expand("resources/database/{ref}/1000G-phase_3.vcf.gz.tbi",ref=config["ref"]["build"])
    params:
        "-p vcf"
    cache: True
    wrapper:
        "master/bio/tabix/index" #no index and 0.78.0

rule tabix_ESP:
    input:
         expand("resources/database/{ref}/ESP6500SI.vcf.gz",ref=config["ref"]["build"])
    output:
         expand("resources/database/{ref}/ESP6500SI.vcf.gz.tbi",ref=config["ref"]["build"])
    params:
        "-p vcf"
    cache: True
    wrapper:
        "master/bio/tabix/index" #no index and 0.78.0     

rule edit_mappability:
    input:
        raw_map=expand("resources/database/{ref}/raw_mappability_{len}.bw",ref=config["ref"]["build"],len=config["kmer"]),
        ref=expand("resources/reference_genome/{ref}/{species}.fasta",ref=config["ref"]["build"],species=config["ref"]["species"]),
        script = "scripts/edit_mappability.R"
    output:
        raw_bed= temp(expand("resources/database/{ref}/raw_mappability_{len}.bed",ref=config["ref"]["build"],len=config["kmer"])),
        bed= temp(expand("resources/database/{ref}/mappability_{len}.bed",ref=config["ref"]["build"],len=config["kmer"])),
        chromsizes= temp(expand("resources/reference_genome/{ref}/sizes.genome",ref=config["ref"]["build"])),
        map= expand("resources/database/{ref}/mappability_{len}.bw",ref=config["ref"]["build"],len=config["kmer"])
    conda:
        "../envs/r4.yaml"
    shell:
        "bigWigToBedGraph {input.raw_map} {output.raw_bed}; "
        "Rscript {input.script} {output.raw_bed} {output.bed}; "
        "faidx {input.ref} -i chromsizes > {output.chromsizes}; "
        "bedGraphToBigWig {output.bed} {output.chromsizes} {output.map}; "
                
##########################
