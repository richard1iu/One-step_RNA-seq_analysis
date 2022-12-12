configfile: "workflow/config/config2.yaml"
data=config["DATA"]
print(data)
#
input_dir="data/paried/alfafa_zm4_reads/"
#
log_dir="log/"
#
fastqc_dir="results/fastqc_results/"
#
fastp_dir="results/trim/fastp/"
trammomatic_dir="results/trim/trimmomatic/"
sickle_dir="results/trim/sickle/"
#  
reference_genome="data/paried/alfafa_zm4_reference_genome/zm-4.genome.fasta"
hisat2_index_dir="data/paried/hisat2_index/"
star_index_dir="data/paried/star_index/"
#
trinity_dir="results/assembly/trinity/"
min_contig_length = 300
max_trinity_memory = "16G"
#
transdecoder_dir = "results/assembly/transdecoder/"
minimum_protein_length = 100
pfam_file="/mnt/f/Database/Pfam/Pfam-A.hmm"

#
bowtie2_dir = "results/assembly/evaluation/bowtie2/"
BUSCO_dir = "results/assembly/evaluation/BUSCO/"
BUSCO_database = "/mnt/f/Database/BUSCO_odb10/embryophyta_odb10"
BUSCO_reads_type = "tran"
#
hisat2_dir="results/align/hisat2/"
star_dir="results/align/star/"
#
reference_gtf = "data/paried/alfafa_zm4_reference_genome/zm-4.gene.gff3"
stringtie_dir = "results/quant/stringtie/"

rule all:
    input:
        expand(stringtie_dir + "quant/{data}/{data}.quant_stringtie.gtf",data = data),
        stringtie_dir+"diff_gene_list.csv"

rule fastqc_fq:
    input:
        r1=input_dir+"{data}_1.fq.gz",
        r2=input_dir+"{data}_2.fq.gz"
    output:
        r1=multiext(fastqc_dir+ "{data}_1",".zip",".html"),
        r2=multiext(fastqc_dir+ "{data}_2",".zip",".html")
    log:
        log_dir+"{data}.fastqc.log"
    params:
        fastqc_dir=fastqc_dir,
        prefix = "{data}"
    threads: 4
    shell:
        "fastqc {input.r1} {input.r2} -o {params.fastqc_dir} 2> {log}"

rule multiqc_fq:
    input:
        fastqc_dir
    output:
        fastqc_dir+"multiqc_report.html"
    threads: 4
    shell:
        "multiqc {input}/. -o {input} "

rule fastp_trim:
    input:
        r1=input_dir+"{data}_1.fq.gz",
        r2=input_dir+"{data}_2.fq.gz"
    output:
        r1_trim=fastp_dir+"{data}_1.trim_fastp.fq.gz",
        r2_trim=fastp_dir+"{data}_2.trim_fastp.fq.gz",
        html=fastp_dir+"{data}.trim_fastp.html",
        json=fastp_dir+"{data}.trim_fastp.json"
    params:
        min_len = 70
    log:
        log_dir+"{data}.trim_fastp.log"
    threads: 4
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
        -o {output.r1_trim} -O {output.r2_trim} \
        --thread {threads} --html {output.html} --json {output.json} 2>{log}
        """

rule trimmomatic_trim:
    input:
        r1=input_dir+"{data}_1.fq.gz",
        r2=input_dir+"{data}_2.fq.gz"
    output:
        r1_trim=trammomatic_dir+"{data}_1.trim_momatic.fq.gz",
        r2_trim=trammomatic_dir+"{data}_2.trim_momatic.fq.gz",
        r1_unpair=trammomatic_dir+"{data}_1.unpair.trim_momatic.fq.gz",
        r2_unpair=trammomatic_dir+"{data}_2.unpair.trim_momatic.fq.gz"
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} {output.r1_trim} {output.r1_unpair} {output.r2_trim} {output.r2_unpair} \
        SLIDINGWINDOW:4:15 MINLEN:25 ILLUMINACLIP:data/adapters/TruSeq3-PE.fa:2:30:10
        """

rule sickle_trim:
    input:
        r1=input_dir+"{data}_1.fq.gz",
        r2=input_dir+"{data}_2.fq.gz"
    output:
        r1_trim=sickle_dir+"{data}_1.trim_sickle.fq.gz",
        r2_trim=sickle_dir+"{data}_2.trim_sickle.fq.gz",
        singles=sickle_dir+"{data}_singles.trim_sickle.fq.gz"
    shell:
        """
        sickle pe -t sanger -f {input.r1} -r {input.r2} \
        -o {output.r1_trim} -p {output.r2_trim} \
        -s {output.singles} \
        -q 35 -l 45
        """

rule trinity_assembly:
    input:
        r1 = fastp_dir+"{data}_1.trim_fastp.fq.gz",
        r2 = fastp_dir+"{data}_2.trim_fastp.fq.gz"
    output:
        outdir = directory(trinity_dir+"{data}_trinity"),
        fasta = trinity_dir+"{data}_trinity.Trinity.fasta",
        stats = trinity_dir+"{data}_trinity.stats"
    params:
        min_length = min_contig_length,
        max_memory = max_trinity_memory,
        label = "{data}"
    threads: 4
    shell:
        """
        Trinity --seqType fq \
        --left {input.r1}  \
        --right {input.r2} \
        --min_contig_length {params.min_length} \
        --CPU {threads} \
        --max_memory {params.max_memory} \
        --output {output.outdir} \
        --full_cleanup

        TrinityStats.pl {output.fasta} > {output.stats}
        """

rule transdecoder_ORFs:
    input:
        trans_fasta = rules.trinity_assembly.output.fasta
    output:
        outdir = directory(transdecoder_dir+"{data}_ORF"),
        pep =  transdecoder_dir+"{data}_ORF/longest_orfs.pep",
        pfam = transdecoder_dir+"{data}_ORF/pfam.domtblout"
    params:
        min_length = minimum_protein_length,
        pfam = pfam_file
    threads: 4
    shell:
        """
        TransDecoder.LongOrfs -m {params.min_length} -t {input} --output_dir {output.outdir}
        hmmscan --cpu {threads} --domtblout {output.pfam} {params.pfam} {output.pep}
        TransDecoder.Predict --cpu {threads} -t {input} --retain_pfam_hits {output.pfam} --output_dir {output.outdir}
        """

rule bowtie2_index:
    input:
        assemblyed_fasta = rules.trinity_assembly.output.fasta
    output:
        index_dir = directory(bowtie2_dir+"{data}"),
        index_file = bowtie2_dir+"{data}_bowtie2_index.txt"
    params:
        prefix="{data}_index"
    threads:
        4
    shell:
        """
        bowtie2-build --threads {threads} {input} {output.index_dir}{params.prefix} 
        touch {output.index_file}
        """

rule bowtie2_align:
    input:
        r1 = fastp_dir+"{data}_1.trim_fastp.fq.gz",
        r2 = fastp_dir+"{data}_2.trim_fastp.fq.gz",
        index_file = rules.bowtie2_index.output.index_dir
    output:
        stats = bowtie2_dir + "{data}_assembly_stats.txt",
        bam   = bowtie2_dir + "{data}.align_bowtie2.bam"
    threads: 4
    shell:
        """
        bowtie2 -p {threads} -q --no-unal -k 20 -x {input.index_file}{data}_index \
        -1 {input.r1} -2 {input.r2}  \
        2> {output.stats} | samtools view -@ {threads} -Sb -o {output.bam}
        """

rule BUSCO_evaluation:
    input:
        assembled_fasta = rules.trinity_assembly.output.fasta
    output:
        outdir = directory(BUSCO_dir+ "{data}")
    threads: 4
    params:
        db = BUSCO_database,
        types = BUSCO_reads_type
    shell:
        """
        busco -i {input} \ # assembled .fasta
        -c {threads} \ # threads
        -o {output} \ # name of output dir 
        -m {param.types} \ # geno/prot/tran
        -l {param.db} 
        --offline 
        """

rule star_index:
    input:
        reference_genome = reference_genome
    output:
        index_dir = directory(star_index_dir),
        index_file = star_index_dir+"star_index.txt"
    threads:
        4
    shell:
        """
        STAR --runMode genomeGenerate --genomeDir {output.index_dir} \
        --genomeFastaFiles {input.reference_genome} \
        --sjdbGTFfile {reference_gtf} --runThreadN {threads}
        touch {output.index_file}
        """

rule star_align:
    input:
        r1 = fastp_dir+"{data}_1.trim_fastp.fq.gz",
        r2 = fastp_dir+"{data}_2.trim_fastp.fq.gz",
        index_file = star_index_dir+"star_index.txt"
    output:
        star_file = star_dir+"{data}/star.txt"
    threads: 4
    shell:
        """
        rm ~/temp -rf
        STAR --runMode alignReads --readFilesCommand zcat \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped None \
        --genomeDir {star_index_dir} \
        --readFilesIn {input.r1} {input.r2} \
        --outFileNamePrefix {star_dir}{data}/{data} \
        --runThreadN {threads} --outTmpDir ~/temp
        touch {output}
        """

rule hisat2_index:
    input:
        reference_genome = reference_genome
    output:
        index_genome = directory(hisat2_index_dir),
        index_file = hisat2_index_dir+"hisat2_index.txt"
    threads:
        4
    shell:
        """
        hisat2-build -p 4 {input} {output.index_genome}hisat2_index
        touch {output.index_file}
        """

rule hisat2_align:
    input:
        r1 = fastp_dir+"{data}_1.trim_fastp.fq.gz",
        r2 = fastp_dir+"{data}_2.trim_fastp.fq.gz",
        index_file = rules.hisat2_index.output.index_file
    output:
        sam = hisat2_dir + "{data}.align_hisat2.sam",
        bam = hisat2_dir + "{data}.align_hisat2.bam"
    threads: 4
    params:
        reference_genome = reference_genome,
    shell:
        """
        hisat2 -p {threads} \
        --dta \
        -x {hisat2_index_dir}hisat2_index \
        -1 {input.r1} -2 {input.r2} \
        -S {output.sam} \
        --dta-cufflinks --no-unal
        samtools sort -@ {threads} -o {output.bam} {output.sam}
        """

rule stringtie_align:
    input:
        bam = rules.hisat2_align.output.bam
    output:
        gtf = stringtie_dir + "align/{data}.align_stringtie.gtf"
    params:
        l = "{data}",
        reference_gtf=reference_gtf
    threads:
        4
    shell:
        """
        stringtie -p {threads} -G {params.reference_gtf} -o {output.gtf} -l {params.l} {input.bam}
        """

rule stringtie_merge:
    input:
        gtf = expand(stringtie_dir + "align/{data}.align_stringtie.gtf",data=data)
    output:
        gtf = "results/quant/stringtie/stringtie_merged.gtf"
    params:
        l = "merge",
        stringtie_dir=stringtie_dir,
        reference_gtf=reference_gtf
    threads:
        4
    shell:
        """
        ls -1 {input.gtf} >> {params.stringtie_dir}/sample_gtf_list.txt
        stringtie --merge -p 4 -G {params.reference_gtf} -o {output.gtf}  {params.stringtie_dir}/sample_gtf_list.txt
        """

rule gffcompare:
    input:
        gtf = rules.stringtie_merge.output.gtf
    output:
        gffcompare = stringtie_dir+"gffcompare/gffcompare.annotated.gtf"
    params:
        stringtie_dir=stringtie_dir,
        reference_gtf=reference_gtf
    threads:
        4
    shell:
        """
        gffcompare -r {params.reference_gtf} -o {params.stringtie_dir}gffcompare/gffcompare {input.gtf}
        """

rule stringtie_quant:
    input:
        bam = rules.hisat2_align.output.bam,
        merged_gtf=rules.stringtie_merge.output.gtf
    output:
        gtf = stringtie_dir + "quant/{data}/{data}.quant_stringtie.gtf"
    threads:
        4
    shell:
        """
        stringtie -e -B -p {threads} -G {input.merged_gtf} -o {output.gtf} {input.bam}
        """

rule stringtie_count:
    input:
        gtf = expand(stringtie_dir + "quant/{data}/{data}.quant_stringtie.gtf",data=data)
    output:
        stringtie_dir+"gene_count_matrix.csv",
        stringtie_dir+"transcript_count_matrix.csv"
    params:
        stringtie_dir=stringtie_dir,
        reference_gtf=reference_gtf
    threads:
        4
    shell:
        """
        ls {input.gtf} > {stringtie_dir}gtf_path.txt
        # find {stringtie_dir} -name *.quant*.gtf > {stringtie_dir}gtf_path.txt
        cat {stringtie_dir}gtf_path.txt | cut -d / -f5 > {stringtie_dir}sampleid.txt
        paste {stringtie_dir}sampleid.txt {stringtie_dir}gtf_path.txt > {stringtie_dir}gtf_list.txt
        rm {stringtie_dir}sampleid.txt {stringtie_dir}gtf_path.txt

        python workflow/scripts/prepDE_py3.py \
        -i {stringtie_dir}gtf_list.txt \
        -t {stringtie_dir}transcript_count_matrix.csv \
        -g {stringtie_dir}gene_count_matrix.csv
        """

rule DEseq2:
    input:
        gene_count = stringtie_dir+"gene_count_matrix.csv",
        transcript_count = stringtie_dir+"transcript_count_matrix.csv"
    output:
        diff_gene_list = stringtie_dir+"diff_gene_list.csv"
    params:
        sample_group = input_dir+"sample_group.txt"
    #conda:
        #"envs/ballgown.yaml"
    script:
        "workflow/scripts/Deseq2.R"