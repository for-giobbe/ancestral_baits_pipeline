#####################################################################################
#                                                                                   #
#                   Script used to design probes for eDNA capture.                  #
#                                                                                   #
#####################################################################################

# to do:
# check numeric flags are integers
# important error: user get it wrong and there are stop codons in the aln
# check that the species in the tree are present in the alignment - PS: species could be added on a backbone

import os
import glob
import shutil
import os.path
import datetime
import argparse
import operator
import subprocess
from os import path
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqIO import FastaIO
from Bio.Align.Applications import MafftCommandline

#################################################################################### parsing arguments and setting variables

def check_positive(value):
    try:
        value = int(value)
        if value <= 0:
            raise argparse.ArgumentTypeError("{} is not a positive integer".format(value))
    except ValueError:
        raise Exception("{} is not an integer".format(value))
    return value

parser = argparse.ArgumentParser(prog='capture', description='design probes for eDNA capture')

parser._optionals.title = "Required Arguments"
optional_args=parser.add_argument_group('Optional Arguments')

parser.add_argument('--out', default="probes", help='basename of oputput file and folders', metavar='')
parser.add_argument('--taxa', help='input taxa', metavar='')
parser.add_argument('--marker', choices=['COI', 'rbcL', 'matK', 'ITS', 'coding' , 'noncoding'], help='marker: can be COI / rbcL / matK / ITS when downloading sequences from BOLD - can be coding/ noncoding when using custom msa', metavar='')
parser.add_argument('--filter_len', '-fl', default=100, help='minimum amminoacids length - defeault is 100', type=check_positive, metavar='')
parser.add_argument('--filter_tax', '-ft', default="", choices=['', 'genus', 'family', 'order'], help='can be genus / family / order - defeault is none', metavar='')
parser.add_argument('--geo', default="", help='geographic location of samples', metavar='')

parser.add_argument('--custom_fas', '-cf', default="", help='custom alignment in fasta format - will skip sequence filtering step', metavar='')
parser.add_argument('--custom_nwk', '-cn', default="", help='custom tree in newick format - will skip phylogenetci inference, a custom alignment needs to be specified', metavar='')

parser.add_argument('--code', required=True, help='genetic code - e.g. 1 for plastid and nuclear, 5 for mitochondrial invertebrate', metavar='')
parser.add_argument('--cores', default=1, help='number of cores used - defeault is 1', metavar='')
parser.add_argument('--verbose', action='store_false', help='keeps temporary folder and files')

args = parser.parse_args()

#################################################################################### errors !!!

if not args.custom_fas and (args.taxa is None):
	print("\n WARNING! When not using custom alignment --taxa flag is required! \n")
	quit()

if not args.custom_fas and (args.marker is None):
	print("\n WARNING! When not using custom alignment --marker flags is required! \n")
	quit()
	
if  args.custom_nwk and (args.custom_fas is None):
	print("\n WARNING! When using a custom tree also a custom alignemtn has to be specified! \n")
	quit()

if path.exists(args.out):
	print("\n WARNING ! An output folder with the same name already exists! \n")
	quit()
else:
	os.makedirs(args.out)

if ( args.custom_fas == "" ) and ( args.custom_nwk != "" ):
	print("\n WARNING ! When specifying a custom phylogeny, also a custo msa has to be specified! \n")
	quit()

os.makedirs(args.out + "/tmp")
os.chdir(args.out + "/tmp")

#################################################################################### define markers

if ( args.marker == "COI" ):
	marker_list=[ "COI" , "COI-5P" , "COI-3P" ]
	coding = True
elif ( args.marker == "rbcL" ):
	marker_list=[ "rbcl" , "RBCL" , "Rbcl" , "rbcla" , "rbcl-a" ]
	coding = True
elif ( args.marker == "matK" ):
	marker_list=[ "matK" , "matk" , "MATK" ]
	coding = True
elif ( args.marker == "ITS" ):
	marker_list=[ "ITS" , "its" ]
	coding = False

if ( args.custom_fas != "" ) and ( args.code != "" ):
	coding = True
else:
	coding = False

#################################################################################### big IF

if (args.custom_fas == ""):

#################################################################################### download from BOLD
	
	if ( args.geo == "" ):
		print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t downloading" , args.marker , "sequences for taxa" , args.taxa)
		subprocess.run(["../../bold-cli" , "-taxon" , args.taxa , "-marker" , args.marker , "-output" , "tmp.bold"] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)
	else:
		print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t downloading" , args.marker , "sequences for taxa" , args.taxa, "from " , args.geo)
		subprocess.run(["../../bold-cli" , "-taxon" , args.taxa , "-marker" , args.marker , "-geo" , args.geo, "-output" , "tmp.bold"] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)
	
	if os.path.getsize("tmp.bold") == 0 :
		print("\n WARNING! The BOLD search for marker" , args.marker , "and taxa" , args.taxa , "returned nothing! \n")
		quit()
	
#################################################################################### filter taxonomy and remove gaps
	
	tmp_1_fasta=[]
	
	if (args.filter_tax != ""):
		print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t filtering sequences with" , args.filter_tax , "identification")
	else:
		print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t no taxonomic filtering")
	
	with open("tmp.bold", 'r', errors='ignore') as file:
		header = file.readline()
		for l in file :
			sl = l.split('\t')
			if (args.filter_tax == "family"):
				if sl[15]:
					header = (">" + sl[0])
					seq = sl[71].replace('-','')
					if sl[69] in marker_list:
						tmp_1_line=[header, seq]
						tmp_1_fasta.append(tmp_1_line)
			elif (args.filter_tax == "genus"):
				if sl[19]:
					header = (">" + sl[0])
					seq = sl[71].replace('-','')
					if sl[69] in marker_list:
						tmp_1_line=[header, seq]
						tmp_1_fasta.append(tmp_1_line)
			elif (args.filter_tax == "species"):
				if sl[21]:
					header = (">" + sl[0])
					seq = sl[71].replace('-','')
					if sl[69] in marker_list:
						tmp_1_line=[header, seq]
						tmp_1_fasta.append(tmp_1_line)
			else:
					header = (">" + sl[0])
					seq = sl[71].replace('-','')
					if sl[69] in marker_list:
						tmp_1_line=[header, seq]
						tmp_1_fasta.append(tmp_1_line)
	
	with open('tmp1.fna', 'w') as tmp_1:
	    for line in tmp_1_fasta:
	        for element in line:
	            tmp_1.write(str(element) + '\n')
	
	if os.path.getsize("tmp1.fna") == 0 :
		print("\n WARNING! No sequence passed the taxonomic filter! \n")
		quit()
	
#################################################################################### collapse identical sequences - can enter custom msa

if (args.custom_fas == ""):
	print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t collapsing identical sequences")
	subprocess.run(["cd-hit", "-i", "tmp1.fna" , "-o" , "tmp2.fna", "-c" , "1.00"] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)
	
if (args.custom_fas != ""):
	custom_msa =  "../../" + args.custom_fas
	print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t custom msa detected")
	print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t collapsing identical sequences")
	subprocess.run(["cd-hit", "-i", custom_msa , "-o" , "tmp2.fna", "-c" , "1.00"] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)

#################################################################################### find orfs

print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t extracting longest orf - minimum of" , args.filter_len)

subprocess.run(["getorf" , "-sequence", "tmp2.fna" , "-outseq" , "tmp1.faa" , "-table" , args.code] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)

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
	if record.id in cds_to_extract and len(record.seq) >= int(args.filter_len):
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

#if num < 4 :
#	print("\n WARNING! Less than 4 sequences passed the length cutoff! \n" , num)
#	quit()

#################################################################################### align aminoacids

print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t aligning sequences")

with open('tmp3.faa', 'w') as tmp3_faa, open('tmp2.faa', 'w') as tmp2_faa:
	os.environ.pop('MAFFT_BINARIES')
	subprocess.call(["mafft" , "--adjustdirection" , "tmp2.faa"], stdout=tmp3_faa, stderr=subprocess.DEVNULL)

#################################################################################### retrotranslate alignment

with open('tmp.aln', 'w') as tmp_aln , open('tmp3.fna', 'w') as tmp3_fna:
	subprocess.run(["pal2nal.pl" , "tmp3.faa" , "tmp3.fna" , "-codontable" , args.code , "-output" , "fasta"], stdout=tmp_aln, stderr=subprocess.DEVNULL)

#################################################################################### fragment size selection

print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t select optimal marker size")

#################################################################################### write partition file and infer tree or detect custom nwk

if (args.custom_nwk != ""):

	custom_tre = "../../" + args.custom_nwk
	print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t custom nwk detected - skipping phylogenetic inference")

	def_tre_file = args.out + ".nwk"
	shutil.copy(custom_tre , def_tre_file)
	
	def_aln_file = args.out + ".aln"
	os.rename('tmp.aln' , def_aln_file)

elif (args.custom_nwk == ""):
	
	print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t inferring tree")
	
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
		subprocess.run(["iqtree" , "-s" , "tmp.aln" , "-nt" , str(args.cores) , "-spp" , "tmp.partitions"], stdout=tmp_nwk, stderr=subprocess.DEVNULL)
		
	def_tre_file = args.out + ".nwk"
	os.rename('tmp.partitions.treefile' , def_tre_file)
	
	def_aln_file = args.out + ".aln"
	os.rename('tmp.aln' , def_aln_file)
	
################################################################################### ancestral sequences

print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t inferring ancestral sequences \n")

icode = int(args.code) - 1

tmp_ctl=[]

aln_line="seqfile = " + def_aln_file
tmp_ctl.append(aln_line)
tmp_ctl.append("outfile = tmp_baseml.out")
tre_line="treefile = " + def_tre_file
tmp_ctl.append(tre_line)
tmp_ctl.append("noisy = 3")
tmp_ctl.append("verbose = 0")
tmp_ctl.append("runmode = 0")
tmp_ctl.append("model = 7")
tmp_ctl.append("Mgene = 0")
tmp_ctl.append("clock = 0")
tmp_ctl.append("fix_kappa = 0")
tmp_ctl.append("kappa = 2.5")
tmp_ctl.append("fix_alpha = 1")
tmp_ctl.append("alpha = 0.")
tmp_ctl.append("Malpha = 0")
tmp_ctl.append("ncatG = 5")
tmp_ctl.append("fix_rho = 1")
tmp_ctl.append("rho = 0.")
tmp_ctl.append("nparK = 0")
tmp_ctl.append("nhomo = 0")
tmp_ctl.append("getSE = 0")
tmp_ctl.append("RateAncestor = 1")
tmp_ctl.append("Small_Diff = 1e-6")
tmp_ctl.append("cleandata = 0")
if coding == True :
	icode_line="icode = " + str(icode)
	tmp_ctl.append(icode_line)
tmp_ctl.append("fix_blength = 2")
tmp_ctl.append("method = 0")

with open('tmp.ctl', 'w') as tmp_ctl_file:
    for line in tmp_ctl:
            tmp_ctl_file.write(line + '\n')
            
with open('tmp.nwk', 'w') as tmp_nwk, open('tmp_baseml.log', 'w') as tmp_baseml_log:
	subprocess.run(["baseml" , "tmp.ctl"], stdout=tmp_baseml_log, stderr=subprocess.DEVNULL)

################################################################################### clean

os.chdir('..')

if args.verbose == False :
	pass
else :
	files = glob.glob("tmp/*")
	for f in files:
		os.remove(f)
	os.rmdir("tmp")
