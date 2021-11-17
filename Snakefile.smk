import os
import sys

##############################################
# Set up config variables
##############################################
# Set up config variables:
if os.path.exists('config.yaml'):
    configfile: 'config.yaml'

fastq_dir = config['fastq_dir']
project_dir = config['project_dir']
ref = config['ref']
st = config['st']
isolate_list_handle = config['isolate_list']
fastq_endings = config['fastq_endings'].split(",")
snippy_minfrac = config['snippy_minfrac']
remove_ref = config['remove_ref']
outgroups = config['outgroups']

iqtree_model = config['model']
iqtree_bb = config['bb']

##############################################
# Control flow for removing ref/outgroups:
##############################################

if remove_ref.split(",")[0] == "" and outgroups.split(",")[0] == "": # only run with ref if neither remove ref nor outgroups set
    remove_seqs_flag = [""]
    remove_seqs_flag_no_dot = [""]
elif not remove_ref.split(",")[0] == "" and outgroups.split(",")[0] == "": # if remove ref set, run without ref (not with)
    remove_seqs_flag = ["no_ref."]
    remove_seqs_flag_no_dot = ["no_ref"]
elif remove_ref.split(",")[0] == "" and not outgroups.split(",")[0] == "":
    remove_seqs_flag = ["with_outgroups.", "without_outgroups."] # if outgroups present, always run with and without them
    remove_seqs_flag_no_dot = ["with_outgroups", "without_outgroups"]
elif not remove_ref.split(",")[0] == "" and not outgroups.split(",")[0] == "":
    remove_seqs_flag = ["no_ref_with_outgroups.", "no_ref_without_outgroups."] # if outgroups present and remove ref, remove ref & run with and without outgroups
    remove_seqs_flag_no_dot = ["no_ref_with_outgroups", "no_ref_without_outgroups"]
else:
    print("Check remove reference/outgroups flags and try again.")
    sys.exit()

##############################################
# Load isolate list and fastq files
##############################################

# Load isolate list:
with open(isolate_list_handle, 'r') as infile1:
    isolate_list = [line.strip() for line in infile1]

# Get list of fastq files:
fastq_files = [os.path.abspath(os.path.join(fastq_dir, fq)) for fq in os.listdir(fastq_dir)]

# Create dict of fastq files for each isolate:
# dict structure: {isolate1: {R1: forward reads, R2: reverse reads}, isolate2: {R1: ..., R2: ...}, ...}
fastq_files_dict = {}
for isolate in isolate_list:
    
    freads = [fastq for fastq in fastq_files if isolate in fastq and fastq_endings[0] in fastq][0]
    rreads = [fastq for fastq in fastq_files if isolate in fastq and fastq_endings[1] in fastq][0]

    fastq_files_dict[isolate] = {'R1':freads, 'R2':rreads}

##############################################
# Set up directories (since snippy-core/iqtree/cfml complain if directory doesn't exist)
##############################################

if not os.path.isdir(os.path.join(project_dir, "snippy-core")):
    os.mkdir(os.path.join(project_dir, "snippy-core"))

if not os.path.isdir(os.path.join(project_dir, "iqtree")): # not sure if iqtree makes an output dir
    os.mkdir(os.path.join(project_dir, "iqtree"))

if not os.path.isdir(os.path.join(project_dir, "cfml")):
    os.mkdir(os.path.join(project_dir, "cfml"))

##############################################
# Define input functions
##############################################

# Def function for getting paired reads:
def get_reads(wildcards):
    return {'R1':fastq_files_dict[wildcards.sample]["R1"], 'R2':fastq_files_dict[wildcards.sample]["R2"]}

##############################################
# All rule
##############################################

# All rule:
rule all:
    input:
        #### snippy: ####
        # expand(os.path.join(project_dir,"raw-snippy-output","{sample}"),sample=isolate_list),
        expand("{output}raw-snippy-output/{sample}", output=project_dir, sample=isolate_list),
        #### snippy-core: ####
        # os.path.join(project_dir,"snippy-core") + "/" + st + ".full.aln",
        expand("{output}snippy-core/{st}.full.aln", output=project_dir, st=st),
        #### snippy-clean_full_aln: ####
        expand("{output}snippy-core/{st}.clean.full.aln", output=project_dir, st=st),
        #### remove_seq: ####
        # expand("{output}snippy-core/{st}.no_ref.clean.full.aln", output=project_dir, st=st),
        # expand("{output}snippy-core/{st}.no_ref_no_outgroups.clean.full.aln", output=project_dir, st=st)
        expand("{output}snippy-core/{st}.{remove_seqs_flag}clean.full.aln", output=project_dir, st=st, remove_seqs_flag=remove_seqs_flag),
        #### iqtree ####
        # expand("{output}iqtree/{st}.{remove_seqs_flag_no_dot}.iqtree", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot)
        expand("{output}iqtree/{st}.{remove_seqs_flag_no_dot}.contree", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot) if len(isolate_list) >= 4 else expand("{output}iqtree/{st}.{remove_seqs_flag_no_dot}.treefile", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot),
        #### cfml ####
        expand("{output}cfml/{st}.{remove_seqs_flag_no_dot}.labelled_tree.newick", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot)

##############################################
# Individual target rules
##############################################

### Run snippy only: ###
rule run_snippy:
    input:
        # expand(os.path.join(project_dir, "raw-snippy-output", "{sample}"), sample=isolate_list),
        expand("{output}raw-snippy-output/{sample}",output=project_dir,sample=isolate_list)

### Run up to snippy core only: ###
rule run_snippy_core:
    input:
        # os.path.join(project_dir,"snippy-core") + "/" + st + ".full.aln"
        expand("{output}snippy-core/{st}.full.aln",output=project_dir,st=st)

### Run up to snippy-clean_full_aln only: ###
rule run_snippy_clean:
    input:
        expand("{output}snippy-core/{st}.clean.full.aln",output=project_dir,st=st)

### Run up to remove ref/outgroups only: ###
rule run_remove_seqs:
    input:
        # expand("{output}snippy-core/{st}.no_ref.clean.full.aln", output=project_dir, st=st),
        # expand("{output}snippy-core/{st}.no_ref_no_outgroups.clean.full.aln", output=project_dir, st=st)
        expand("{output}snippy-core/{st}.{remove_seqs_flag}clean.full.aln", output=project_dir, st=st, remove_seqs_flag=remove_seqs_flag),

### Run up to iqtree only: ###
rule run_iqtree:
    input:
        # expand("{output}iqtree/{st}.{remove_seqs_flag_no_dot}iqtree", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot)
        expand("{output}iqtree/{st}.{remove_seqs_flag_no_dot}.contree",output=project_dir,st=st,remove_seqs_flag_no_dot=remove_seqs_flag_no_dot) if len(isolate_list) >= 4 else expand("{output}iqtree/{st}.{remove_seqs_flag_no_dot}.treefile",output=project_dir,st=st,remove_seqs_flag_no_dot=remove_seqs_flag_no_dot)

### Run up to cfml only: ###
rule run_cfml:
    input:
        expand("{output}cfml/{st}.{remove_seqs_flag_no_dot}.labelled_tree.newick",output=project_dir,st=st,remove_seqs_flag_no_dot=remove_seqs_flag_no_dot)

##############################################
# Individual rules
##############################################

### snippy rule: ###
rule snippy:
    input:
        # R1 = "/home/conrad/test-snakemake/fastq-files/{sample}_1.fastq.gz",
        # R2 = "/home/conrad/test-snakemake/fastq-files/{sample}_2.fastq.gz"
        unpack(get_reads)
    wildcard_constraints: # I think the issue here is that maybe because in input, fastq-files/ prepends {sample}, while it does not in output (i.e. input and output are in different directories), which may confuse snakemake and make it try to set the wildcard equal to the difference between the two directories (i.e. add fastq-files/ to isolate name).
        sample = '[a-zA-Z0-9_-]+' # doesn't work without this...expands sample to include fastq-files/ for some reason...no idea why when sample can only take on values in isolate_list in 'rule snippy' above...
    params:
        reference=ref,
        minfrac=snippy_minfrac
    output:
        # directory(os.path.join(project_dir, "raw-snippy-output", "{sample}"))
        directory("{output}raw-snippy-output/{sample}")
    shell:
        "snippy --R1 {input.R1} --R2 {input.R2} --ref {params.reference} --outdir {output} --cpus 8 --minfrac {params.minfrac}"

### snippy-core rule: ###
rule snippy_core:
    input:
        expand("{output}raw-snippy-output/{sample}", output=project_dir, sample=isolate_list)
    output:
        "{output}snippy-core/{st}.full.aln"
    params:
        reference=ref,
        output_prefix="{output}/snippy-core/{st}"
    wildcard_constraints:
        # needed otherwise snakemake confuses the output here with the output (no ref output) from
        # the remove ref step...idk why, because the filenames are different (one has .no_ref, the
        # other doesn't, but apparently the wildcards get confused...
        st = '[a-zA-Z0-9_]+'
    shell:
        "snippy-core --prefix {params.output_prefix} --ref {params.reference} {input}"

### snippy clean rule: ###
rule snippy_clean:
    input:
        "{output}snippy-core/{st}.full.aln"
    output:
        "{output}snippy-core/{st}.clean.full.aln"
    shell:
        "snippy-clean_full_aln {input} > {output}"

### remove seqs rule: ###
rule remove_seqs:
    input:
        expand("{output}snippy-core/{st}.clean.full.aln", output=project_dir, st=st) # Output here is not the output of any other script, so wildcards need to be defined here.
    output:
        expand("{output}snippy-core/{st}.{remove_seqs_flag}clean.full.aln", output=project_dir, st=st, remove_seqs_flag=remove_seqs_flag)
    params:
        remove_ref = remove_ref,
        remove_outgroups = outgroups
    script:
        "/home/conrad/test-snakemake/remove_seqs_function.py"

### iqtree rule: ###
iqtree_params_string = ""
if len(isolate_list) >= 4:
    iqtree_params_string = f"-m {iqtree_model} -bb {iqtree_bb} -nt AUTO"
elif len(isolate_list) < 4:
    iqtree_params_string = f"-m {iqtree_model} -nt AUTO"

rule iqtree:
    input:
        "{output}snippy-core/{st}.{remove_seqs_flag_no_dot}.clean.full.aln"
    output:
        "{output}iqtree/{st}.{remove_seqs_flag_no_dot}.contree" if len(isolate_list) >= 4 else "{output}iqtree/{st}.{remove_seqs_flag_no_dot}.treefile"
    params:
        iqtree_params = iqtree_params_string,
        iqtree_prefix = "{output}iqtree/{st}.{remove_seqs_flag_no_dot}"
    shell:
        "iqtree -s {input[0]} {params.iqtree_params} -pre {params.iqtree_prefix}"

### ClonalFrameML rule: ###
# cfml_ending = ""
# if len(isolate_list) >= 4:
#     cfml_ending = ".contree"
# elif len(isolate_list) < 4:
#     cfml_ending = ".treefile"
#
rule clonalframeml:
    input:
        tree = "{output}iqtree/{st}.{remove_seqs_flag_no_dot}.contree" if len(isolate_list) >= 4 else "{output}/iqtree/{st}.{remove_seqs_flag_no_dot}.iqtree",
        fasta = "{output}snippy-core/{st}.{remove_seqs_flag_no_dot}.clean.full.aln"
    output:
        "{output}cfml/{st}.{remove_seqs_flag_no_dot}.labelled_tree.newick"
    params:
        cfml_prefix = "{output}cfml/{st}.{remove_seqs_flag_no_dot}"
    shell:
        "ClonalFrameML {input.tree} {input.fasta} {params.cfml_prefix} -em true -show_progress true"