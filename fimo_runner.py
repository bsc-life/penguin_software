import os
import sys
import StringIO
from Bio import SeqIO
import subprocess
import IPython
import shutil



####################
#### please make sure you have installed FIMO and included it your PATH
#### To run:
#### python fimo_runner.py sequences.fasta
####################
dnaSeq = {}
in_arg = sys.argv[1]
input_handle = open(in_arg, "rU")
for record in SeqIO.parse(input_handle, "fasta"):
    if record.id not in dnaSeq:
        dnaSeq[record.id] = record

out = os.path.join("src/EPIN_reconstruction_data/", "FIMO_outputs_{}".format(in_arg.replace(".fasta", "").replace(".txt", "")))
if not os.path.exists(out):
    os.makedirs(out)

motifs = os.path.join("src/EPIN_reconstruction_data/", "JASPAR_motifs_2020_human_redundant_meme.txt")


for record , seq in dnaSeq.iteritems():

    output_handle = open("{}_dna_tmp2.fasta".format(in_arg), "w")
    SeqIO.write(seq, output_handle, "fasta")
    output_handle.close()

    
    command = "fimo --verbosity 1 -norc -oc {}_fimo_out2 {} {}_dna_tmp2.fasta".format(in_arg, motifs,  in_arg)
    
    p = subprocess.Popen(command, shell=True)
    p.communicate()

    
    p = subprocess.Popen("mv {}_fimo_out2/fimo.tsv {}/{}.tsv".format(in_arg, out, record.replace(":","_")), shell=True)
    p.communicate()
    os.remove("{}_dna_tmp2.fasta".format(in_arg))
    shutil.rmtree("{}_fimo_out2".format(in_arg))
