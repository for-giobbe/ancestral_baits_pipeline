#####################################################################################
#                                                                                   #
#               This script is used to design baits for eDNA capture.               #
#                                                                                   #
##################################################################################### 

import os
import glob
import argparse
import operator
import subprocess
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqIO import FastaIO
from Bio.Align.Applications import MafftCommandline

#################################################################################### parsing arguments and setting variables

parser = argparse.ArgumentParser(description='design baits for eDNA capture')

parser.add_argument('--input_taxa', required=True, help='Input taxa.')
parser.add_argument('--marker', required=True, help='marker (e.g. COI-5P)')
parser.add_argument('--gen_code', required=True, help='Genetic code (e.g. 1 for plastid and nuclear, 5 for mitochondrial invertebrate)')
parser.add_argument('--basename', nargs='?', default="output", help='base name of the oputputs')
parser.add_argument('--min_length', nargs='?', default=100, help='Minimum amminoacids length - defeault is 200')
parser.add_argument('--taxonomy_filter', nargs='?', default="", help='can be genus / family / order - defeault is none')
parser.add_argument('--threads', nargs='?', default=1, help='number of threads used - defeault is 1')
parser.add_argument('--location', nargs='?', default="", help='geographic location of samples')

args = parser.parse_args()

#################################################################################### download from BOLD

if ( args.location == "" ):
	subprocess.run(["./bold-cli" , "-taxon" , args.input_taxa , "-output" , "tmp.bold"] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)
	print("\n\t downloading sequences for taxa" , args.input_taxa)
else:
	subprocess.run(["./bold-cli", "-taxon", args.input_taxa, "-geo" , args.location, "-output" , "tmp.bold"] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)
	print("\n\t downloading sequences for taxa" , args.input_taxa, "from location:" , args.location)

#################################################################################### filter taxonomy and remove gaps

tmp_1_fasta=[]

if (args.taxonomy_filter != ""):
	print("\n\t filtering sequences with" , args.taxonomy_filter , "identification")
else:
	print("\n\t no taxonomic filtering")

with open("tmp.bold", 'r') as file:
	header = file.readline()
	for l in file :
		sl = l.split('\t')
		if (args.taxonomy_filter == "family"):
			if sl[15]:
				header = (">" + sl[0])
				seq = sl[71].replace('-','')
				tmp_1_line=[header, seq]
				tmp_1_fasta.append(tmp_1_line)
		elif (args.taxonomy_filter == "genus"):
			if sl[19]:
				header = (">" + sl[0])
				seq = sl[71].replace('-','')
				tmp_1_line=[header, seq]
				tmp_1_fasta.append(tmp_1_line)
		elif (args.taxonomy_filter == "species"):
			if sl[21]:
				header = (">" + sl[0])
				seq = sl[71].replace('-','')
				tmp_1_line=[header, seq]
				tmp_1_fasta.append(tmp_1_line)
		else:
				header = (">" + sl[0])
				seq = sl[71].replace('-','')
				tmp_1_line=[header, seq]
				tmp_1_fasta.append(tmp_1_line)
                
with open('tmp1.fna', 'w') as tmp_1:
    for line in tmp_1_fasta:
        for element in line:
            tmp_1.write(str(element) + '\n')

#################################################################################### collapse identical sequences

print("\n\t collapsing identical sequences")

subprocess.run(["cd-hit", "-i", "tmp1.fna" , "-o" , "tmp2.fna", "-c" , "1.00"] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)

#################################################################################### find orfs

print("\n\t extracting longest orf - minimum of" , args.min_length)

subprocess.run(["getorf" , "-sequence", "tmp2.fna" , "-outseq" , "tmp1.faa" , "-table" , args.gen_code] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)

#################################################################################### find longest orfs and length cutoff

tmp1_faa=SeqIO.parse("tmp1.faa","fasta")
tmp2_faa=open("tmp2.faa",'w')
tmp3_fna=open("tmp3.fna",'w')

uniq = []
cds_to_extract = []

for record in tmp1_faa:
    sample = (record.id.split('_')[0])
    cds_id = (record.id.split('_')[1])
    length = (len(record.seq))
    if not sample in uniq:
        uniq.append(sample)

for sample in uniq:
    cds = {}
    for record in SeqIO.parse("tmp1.faa", "fasta"):
        if (record.id.split('_')[0]) == sample:
                cds[(record.id.split('_')[1]  + " ")] = len(record.seq)
    longest_cds=max(cds, key=cds.get)
    name=((sample+"_"+longest_cds).replace(' ',''))
    cds_to_extract.append(name)

selected_seqs_aa = list()
selected_seqs_header = list()

for record in SeqIO.parse("tmp1.faa", "fasta"):
	if record.id in cds_to_extract and len(record.seq) >= int(args.min_length):
		selected_seqs_header.append(record.id.split('_')[0])     
		record.id = record.id.split('_')[0]
		record.description = ""
		selected_seqs_aa.append(record)
SeqIO.write(selected_seqs_aa , tmp2_faa , "fasta")

selected_seqs_nt = list()

for record in SeqIO.parse("tmp2.fna", "fasta"):
	if record.id in selected_seqs_header:
		selected_seqs_nt.append(record)
SeqIO.write(selected_seqs_nt , tmp3_fna , "fasta")   

#################################################################################### align aminoacids

print("\n\t aligning aminoacids")

with open('tmp3.faa', 'w') as tmp3_faa, open('tmp2.faa', 'w') as tmp2_faa:
	subprocess.call(["mafft" , "--adjustdirection" , "tmp2.faa"], stdout=tmp3_faa, stderr=subprocess.DEVNULL)

#################################################################################### retrotranslate alignment

with open('tmp.aln', 'w') as tmp_aln , open('tmp3.fna', 'w') as tmp3_fna:
	subprocess.run(["pal2nal.pl" , "tmp3.faa" , "tmp3.fna" , "-codontable" , args.gen_code , "-output" , "fasta"], stdout=tmp_aln, stderr=subprocess.DEVNULL)

#################################################################################### write partition file and infer tree

print("\n\t inferring tree \n")

tmp_partitions=[]
filename = "tmp.aln"
format = "fasta"
tmp_aln = AlignIO.read(filename, format)
tmp_partitions.append("DNA, st = 1-%i\\3" % tmp_aln.get_alignment_length())
tmp_partitions.append("DNA, nd = 2-%i\\3" % tmp_aln.get_alignment_length())
tmp_partitions.append("DNA, rd = 3-%i\\3" % tmp_aln.get_alignment_length())

with open('tmp.partitions', 'w') as tmp_partitions_file:
    for line in tmp_partitions:
            tmp_partitions_file.write(line + '\n')

with open('tmp.nwk', 'w') as tmp_nwk :
	subprocess.run(["iqtree" , "-s" , "tmp.aln" , "-nt" , str(args.threads) , "-spp" , "tmp.partitions"], stdout=tmp_nwk, stderr=subprocess.DEVNULL)

################################################################################### clean

def_tre_file = args.basename + ".nwk"
os.rename('tmp.partitions.treefile' , def_tre_file)

def_aln_file = args.basename + ".aln"
os.rename('tmp.aln' , def_aln_file)

tmp_files = glob.glob('tmp*.*')

for filePath in tmp_files:
	os.remove(filePath)