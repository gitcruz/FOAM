from datetime import datetime
import re
import os
import subprocess # important to call pre-existing modules

#1- allows to call subworkflows (modules)
module polish_workflow:
  snakefile: "../modules/polish.rules.smk"
#evaluate_assemblies.rules.smk
module evaluate_assemblies_workflow:
  snakefile: "../modules/evaluate_assemblies.rules.smk"


#2- main snakemake rules
date = datetime.now().strftime('%Y%m%d.%H%M%S')
#get the bin path for the snakemake file
bin_path = sys.path[0]
#make scripts and utils relative to bin
scripts_dir = bin_path + "/../scripts/"
utils_dir = bin_path + "/../utils/"
keepfiles = False
work_dir = os.getcwd() + "/"
logs_dir = work_dir + "logs/"
if not os.path.exists(logs_dir):
 os.makedirs(logs_dir)

rule format_mito_ref:
  input:
    reference = config["data"]["reference_fasta"],
  output:
    mitoref = work_dir + "mitogenome_ref/MT_REF.scaffolds.fa",
  params:
    scripts = scripts_dir,
  log:
    logs_dir+ str(date) + ".j%j.format_mito_ref.out",
    logs_dir+ str(date) + ".j%j.format_mito_ref.err",
  threads: 1
  shell:
    # Just rename it and format the fasta as multiline
    "mkdir -p mitogenome_ref; "
    "cd mitogenome_ref; "
    "ln -s {input.reference} . ;"
    "{params.scripts}rename_fasta_seq.pl -f  {input.reference} -n MT_REF | {params.scripts}FastaToTbl | {params.scripts}TblToFasta > {output.mitoref} ;"
    "echo 'reference mitogenome has been formatted into a multiline fasta'; "

rule align_ont_to_mitoref:
  input:
    genome = rules.format_mito_ref.output.mitoref,
    reads = config["data"]["ont"]
  output:
    ont_mito_fq = work_dir + "mappings2ref/ont_mito.fastq"
  params:
    align_opts = config["parameters"]["minimap2_opts"],
    tmp = work_dir + "tmp/minimap2.sam",
    scripts = scripts_dir,
    min_match = config["parameters"]["ontfilter_minmatch"],
    min_qual = config["parameters"]["ontfilter_minq"],
  conda:
    "../envs/fa5f7adcb4c59485a3a283d809e1402e.yaml"
  threads: 8
  log:
    logs_dir+ str(date) + ".j%j.align_ont_to_mito_ref.out",
    logs_dir+ str(date) + ".j%j.align_ont_to_mito_ref.err",
  shell:
    "mkdir -p tmp; "
    "mkdir -p mappings2ref; "
    "cd  mappings2ref; "
    #getting REFNAME for fishing reads
    "REFNAME=$( grep \'>\' {input.genome} | sed \'s/>//g\'); "
    "echo \"mitochondrial reference name is $REFNAME\";"

    # pipe to get a sam with the aligned reads to the mitogenome reference calling the md tag
    "minimap2 {params.align_opts} -t {threads} {input.genome} {input.reads} | samtools calmd -S - {input.genome} > {params.tmp};" # specify the output sam is now in the tmpdir

    # CAPTURE READS 1. Using defaults here. my $minmatch = 800 (hard-coded in perl-script) and $minqual = 12 (parameter in perlscript) 
    "perl --version ;"
    "ls -l {params.scripts}sam_id_filter_v02.pl ;"
    "{params.scripts}sam_id_filter_v02.pl -ref $REFNAME -m {params.min_match} -q {params.min_qual} <  {params.tmp} > ont.clean.fastq 2> {output.ont_mito_fq} ;"
    
rule flye_assembly_meta:
  input:
    reads = rules.align_ont_to_mitoref.output.ont_mito_fq
  output:
    assembly = work_dir+"flye_meta/flye_meta.fa",
    info = work_dir+"flye_meta/flye_meta_info.txt",
    png = work_dir+"flye_meta/flye_meta.png",
    circular_contig = work_dir+"flye_meta/circular_contig.fa",
    circular_bed = work_dir+"flye_meta/circular_contig.bed",
  params:
    scripts = scripts_dir,
    out_dir = work_dir+"flye_meta/",
    assembly_opts = config["parameters"]["flye_opts"],
    ref_len = config["parameters"]["reference_length"],
    len_var = config["parameters"]["mito_length_variation"],
  conda:
    "../envs/36b7263cd61f32e4e9f7db777993cd52.yaml"   
  threads: 16
  log:
    logs_dir+ str(date) + ".j%j.flye_meta.out",
    logs_dir+ str(date) + ".j%j.flye_meta.err",
  shell:
    "mkdir -p {params.out_dir}out; "
    "cd {params.out_dir}; "
    #added --resume options to avoid redoing flye if computed successfuly (see assembly_opts under rule params)
    "echo 'Running command: flye {params.assembly_opts} -g {params.ref_len} -o {params.out_dir}out -t {threads} --nano-raw {input.reads}'; "
    "flye {params.assembly_opts} -g {params.ref_len} -o {params.out_dir}out -t {threads} --nano-raw {input.reads}; "
    "ln -s {params.out_dir}out/assembly.fasta {output.assembly}; "
    "ln -s {params.out_dir}out/assembly_info.txt {output.info}; "
    "export XDG_RUNTIME_DIR=$PWD ;" # to avoid Bandage execution error
    "Bandage image {params.out_dir}out/assembly_graph.gfa {output.png}; "    
    "echo 'Selecting the circular contig with length between {params.ref_len}-{params.len_var} and {params.ref_len}+{params.len_var}'; "
    "circular=$(cat {output.info} | gawk \'{{ if ( ($4 == \"Y\") &&  ($2 >= ({params.ref_len}-{params.len_var})) && ($2 <= ({params.ref_len}+{params.len_var})) )   print $1}}\') ; " #double brackets to make it work
    "echo circular contig\(s\)\: $circular ;"
    "{params.scripts}select_contig.pl -f {output.assembly} -t $circular > {output.circular_contig} ; "
    # want to store the coordinates and name of coordinate contig to pass it to samtools in next rule align_illumina_flye_meta
    "{params.scripts}fastalength {output.circular_contig} | gawk \'{{ print $2\"\t\"0\"\t\"$1}}\'  > {output.circular_bed} ; "

use rule align_illumina from evaluate_assemblies_workflow with:
  input: 
    genome = rules.flye_assembly_meta.output.assembly,
    reads = [config["data"]["ill_1"],config["data"]["ill_2"]],
  output:
    mapping = work_dir + "illumina2flye/illumina_mito.bam",
    stats = work_dir + "illumina2flye/illumina_mito.stats.txt",
  params:
    options = "",
  log:
    logs_dir+ str(date) + ".j%j.align_illumina2flye.out",
    logs_dir+ str(date) + ".j%j.align_illumina2flye.err",  
  threads: 16

rule select_illumina_from_circular:
  input:
    genome = rules.flye_assembly_meta.output.circular_contig,
    bam = rules.align_illumina.output.mapping,
    bed=rules.flye_assembly_meta.output.circular_bed,
    circular=rules.flye_assembly_meta.output.circular_contig,  
  output:
    bam = work_dir + "illumina2flye/illumina_mito.circular.nodups.bam",
    stats = work_dir + "illumina2flye/illumina_mito.circular.nodups.stats",
    fastq1 = work_dir + "illumina2flye/illumina_mito.circular.nodups.1.fastq",
    fastq2 = work_dir + "illumina2flye/illumina_mito.circular.nodups.2.fastq",
  params:
    options ="-F 1024 -f 3",  
  conda:
    "../envs/17e4d03683d3b1ea1e781fcb2cfef059.yaml"
  threads: 4
  log:
    logs_dir+ str(date) + ".j%j.select_illumina_from_circular.out",
    logs_dir+ str(date) + ".j%j.select_illumina_from_circular.err",
  shell:
    "cd  illumina2flye; "
    #1.params options simultaneously select to exclude PCR/optical dups & keep properly paired alignments
    "samtools view -b -t {threads} {params.options} -L {input.bed} {input.bam} > {output.bam}; "
    "samtools stats {output.bam} > {output.stats}; "
    #2.prepare for igv visualization:
    #index the bam file 
    "samtools index {output.bam} ; "
    #fai index of circular ref
    "ln -s {input.circular} . ; " 
    "samtools faidx circular_contig.fa ; "
    #3. pull out fastq with illumina reads (mitochondrial)
    "samtools fastq -1 {output.fastq1} -2 {output.fastq2} {output.bam} ; " 

#new rule to nexpolish the circular contig using illumina fastqs
rule nextpolish_circular_contig:
  input:
    genome = rules.flye_assembly_meta.output.circular_contig,
    fastq1 = rules.select_illumina_from_circular.output.fastq1,
    fastq2 = rules.select_illumina_from_circular.output.fastq2,
  output:
    polished = work_dir + "nextpolish/out-np-ill2/genome.nextpolish.fasta",
  params:
    scripts = scripts_dir,
  threads: 16
  log:
    logs_dir+ str(date) + ".j%j.nextpolish_circular_contig.out",
    logs_dir+ str(date) + ".j%j.nextpolish_circular_contig.err",
  shell:
    "cd  nextpolish ; "
   
    #1. Load NextPolish module
    "export MODULEPATH=/software/assembly/easybuild/modules/all/:$MODULEPATH ; "
    "module load  NextPolish/1.4.1-GCC-11.2.0 ; "
    "nextPolish --version ; "
    "echo 'Polishing {input.genome} with 2 rounds of illumina short reads' ; "
  
    #2. Link target genome
    "ln -s  {input.genome} . ; "

    #3. create or link file of file names illumina
     "echo {input.fastq1} > illumina.fofn ; "
     "echo {input.fastq2} >> illumina.fofn ; "
   
    #4. copy config here because the locations are relative to the config file location 
     "cp {params.scripts}mito-nextpolish-ill2.cfg . ; "
   
    #5. Run nextPolish:
    "nextPolish mito-nextpolish-ill2.cfg ; "

# new rule to orient the polsihed mitogenome based on the mitos annotation. uses script orient_mitogenome_v1.pl 
rule orient_polished_contig:
  input:
    genome = rules.nextpolish_circular_contig.output.polished,
  output:
    oriented = work_dir + "orient_mitogenome/out/circular_polished.oriented.fa",
    annotation = work_dir + "orient_mitogenome/out/annotation/result.gff",
    formatted = work_dir + "orient_mitogenome/out/formatted.fa",
  params:
    scripts = scripts_dir,
    refseq_dir = utils_dir + "refseq_dir/",
    genetic_code = config["annotation"]["genetic_code"],
    mitos_options = config["annotation"]["mitos_options"],
    refseq_db = config["annotation"]["refseq_database"],
    tolid = config["parameters"]["tolid"],
  conda:
    "../envs/mitos.yaml"
  threads: 8
  log:
    logs_dir+ str(date) + ".j%j.orient_mitogenome.out",
    logs_dir+ str(date) + ".j%j.orient_mitogenome.err",
  shell:
    "cd  orient_mitogenome ; "
   
    #1. link polished contig 
    "ln -s {input.genome} circular_polished.fa ; "

    #2. Annotation 1
    "echo 'testing mitos conda environment' ; "
    "runmitos.py --version ; "
    "mkdir -p annotation_circular ; "
    "runmitos.py -i circular_polished.fa -c {params.genetic_code} {params.mitos_options} -r {params.refseq_db} -R {params.refseq_dir} -o annotation_circular ; "
    
    #3. Orient based in Annotation set trnF start
    "perl {params.scripts}orient_mitogenome_v1.pl  -f circular_polished.fa -b annotation_circular/result.bed ; "
    "mkdir -p out/annotation ; "
    "mv circular_polished.oriented.fa out/; "
    #4. rename and format to ensure multiline fasta
     "cd out/; "
     "{params.scripts}rename_fasta_seq.pl -f  {output.oriented} -n {params.tolid}_MT | {params.scripts}FastaToTbl | {params.scripts}TblToFasta > {output.formatted} ;"
    # 5. Reannotate formatted mitogenome
      "runmitos.py -i {output.formatted} -c {params.genetic_code} {params.mitos_options} -r {params.refseq_db} -R {params.refseq_dir} -o annotation ; "
    # 6. format annotation gff and copy final fasta 
      "cd ../ ; "
      "cp {output.formatted}  {params.tolid}_MT.fa ; "
      "{params.scripts}parse_mitos_gff.pl -g {output.annotation} -p MT{params.tolid} > MT{params.tolid}.gff3 ; "

rule orient_reference_MT:
  input:
    genome = rules.format_mito_ref.output.mitoref,
  output:
    oriented = work_dir + "mitogenome_ref/out/MT_REF.scaffolds.oriented.fa", 
    annotation = work_dir + "mitogenome_ref/out/annotation_2/result.gff",
    formatted = work_dir + "mitogenome_ref/out/MT_REF.formatted.fa",
  params:
    scripts = scripts_dir,
    refseq_dir = utils_dir + "refseq_dir/",
    genetic_code = config["annotation"]["genetic_code"],
    mitos_options = config["annotation"]["mitos_options"],
    refseq_db = config["annotation"]["refseq_database"],
    tolid = config["parameters"]["tolid"],
  conda:
    "../envs/mitos.yaml"
  threads: 8
  log:
    logs_dir+ str(date) + ".j%j.orient_reference_MT.out",
    logs_dir+ str(date) + ".j%j.orient_reference_MT.err",
  shell:
    "cd  mitogenome_ref ; "
   
    # Do not link polished contig 
    
    #1. Annotation 1
    "echo 'testing mitos conda environment' ; "
    "runmitos.py --version ; "
    "mkdir -p annotation_1; "
    "runmitos.py -i {input.genome} -c {params.genetic_code} {params.mitos_options} -r {params.refseq_db} -R {params.refseq_dir} -o annotation_1 ; "
    
    #3. Orient based in Annotation set trnF start
    "perl {params.scripts}orient_mitogenome_v1.pl  -f {input.genome} -b annotation_1/result.bed ; "
    "mkdir -p out/annotation_2 ; "
    "mv MT_REF.scaffolds.oriented.fa out/; "
    #4. rename and format to ensure multiline fasta
     "cd out/; "
     "{params.scripts}rename_fasta_seq.pl -f  {output.oriented} -n {params.tolid}_MT | {params.scripts}FastaToTbl | {params.scripts}TblToFasta > {output.formatted} ;"
    # 5. Reannotate formatted mitogenome
      "runmitos.py -i {output.formatted} -c {params.genetic_code} {params.mitos_options} -r {params.refseq_db} -R {params.refseq_dir} -o annotation_2 ; "
    # 6. format annotation gff and copy final fasta 
      "cd ../ ; "
      "cp {output.formatted}  {params.tolid}_MT.fa ; "
      "{params.scripts}parse_mitos_gff.pl -g {output.annotation} -p MT{params.tolid} > MT{params.tolid}.gff3 ; "