from Bio import SeqIO
import sys

def remove_seqs(seq_to_remove, fasta, output_fasta):
    records = list(SeqIO.parse(fasta, "fasta"))
    records2 = [i for i in records if i.id not in seq_to_remove.split(",")]

    with open(output_fasta, "w") as outfile:
        SeqIO.write(records2, outfile, "fasta")

# remove ref only:
# remove_seqs(snakemake.params[0], snakemake.input[0], snakemake.output[0])
#
# # remove ref and outgroups:
# remove_seqs(snakemake.params[0] + "," + snakemake.params[1], snakemake.input[0], snakemake.output[1])

# Control flow for removing ref/outgroups:
if snakemake.params[0].split(",")[0] == "" and snakemake.params[1].split(",")[0] == "": # only run with ref if neither remove ref nor outgroups set
    print("Not removing anything, running with clean.{st}.full.aln...")
    remove_seqs_flag = [""]
    remove_seqs("", snakemake.input[0], snakemake.output[0])

elif not snakemake.params[0].split(",")[0] == "" and snakemake.params[1].split(",")[0] == "": # if remove ref set, run without ref (not with)
    print("Removing ref and running without ref only...")
    remove_seqs_flag = ["no_ref."]
    remove_seqs(snakemake.params[0], snakemake.input[0], snakemake.output[0])

elif snakemake.params[0].split(",")[0] == "" and not snakemake.params[1].split(",")[0] == "":
    print("Keeping ref and running with/without outgroups...")
    remove_seqs_flag = ["with_outgroups.", "no_outgroups."] # if outgroups present, always run with and without them
    remove_seqs("", snakemake.input[0], snakemake.output[0])
    remove_seqs(snakemake.params[1], snakemake.input[0], snakemake.output[1])

elif not snakemake.params[0].split(",")[0] == "" and not snakemake.params[1].split(",")[0] == "":
    print("Removing ref and running with/without outgroups...")
    remove_seqs_flag = ["no_ref_with_outgroups.", "no_ref_no_outgroups."] # if outgroups present and remove ref, remove ref & run with and without outgroups
    remove_seqs(snakemake.params[0], snakemake.input[0], snakemake.output[0])
    remove_seqs(snakemake.params[0] + "," + snakemake.params[1], snakemake.input[0], snakemake.output[1])

else:
    print("Check remove reference/outgroups flags and try again.")
    sys.exit()
