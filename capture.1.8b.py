#####################################################################################
#                                                                                   #
#                   Script used to design probes for eDNA capture.                  #
#                                                                                   #
#####################################################################################

# to do:
# branchlength in constrain tree
# stop codon warning
# less than 4 sequences warning
# print non-coding status in stdout messages

import os
import re
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

parser.add_argument('--out', default="probes", help='basename of oputput file and folders')
parser.add_argument('--taxa', help='input taxa')
parser.add_argument('--marker', required = True, choices=['COI-5P', 'COI-3P', 'rbcL', 'matK', 'ITS', 'coding' , 'noncoding'], help='marker: can be COI / rbcL / matK / ITS when downloading sequences from BOLD - can be coding/ noncoding when using custom fasta file')
parser.add_argument('--filter_len', '-fl', default=100, help='minimum amminoacids / nucleotides length - defeault is 100', type=check_positive)
parser.add_argument('--filter_tax', '-ft', default="", choices=['genus', 'family', 'order'], help='can be genus / family / order - defeault is none')
parser.add_argument('--geo', default="", help='geographic location of samples')
parser.add_argument('--custom_fas', '-cf', default="", help='custom alignment in fasta format')
parser.add_argument('--custom_nwk', '-cn', default="", help='custom tree in newick format - a custom alignment needs to be specified: skips phylogenetic inference if there is 1:1 match betwenn tree and msa, otherwise will use it as a constrain')
parser.add_argument('--code', help='genetic code - e.g. 1 for plastid and nuclear, 5 for mitochondrial invertebrate')
parser.add_argument('--cores', type=check_positive, default=1, help='number of cores used - defeault is 1')
parser.add_argument('--verbose', action='store_false', help='keeps temporary folder and files')
parser.add_argument('--erase', action='store_true', help='erases and rewrites an output folder with same name')
parser.add_argument('--distance',  default=0.01, help='distance from nodes - max number of substitutions per site')

args=parser.parse_args()

#################################################################################### errors !!!

if not args.custom_fas and (args.taxa is None):
	print("\n WARNING! Either a custom alignment (--custom_fas) or a lineage do download from BOLD (--taxa) is required! \n")
	quit()

if not args.custom_fas and (args.marker is None):
	print("\n WARNING! When not using custom alignment --marker flags is required! \n")
	quit()
	
if  args.custom_nwk and (args.custom_fas is None):
	print("\n WARNING! When using a custom tree also a custom alignment has to be specified! \n")
	quit()

if args.erase == False and path.exists(args.out):
	print("\n WARNING! An output folder with the same name already exists! \n")
	quit()
elif args.erase == True and path.exists(args.out):
	shutil.rmtree(args.out)
	os.makedirs(args.out)
else:
	os.makedirs(args.out)

if ( args.custom_fas == "" ) and ( args.custom_nwk != "" ):
	print("\n WARNING! When specifying a custom phylogeny, also a custom msa has to be specified! \n")
	quit()
	
if (args.code != None):
	if int(args.code) not in range (1, 33):
		print("\n WARNING! Unknown geneitc code has been specified - please refer to NCBI! \n")
		quit()

os.makedirs(args.out + "/tmp")
os.chdir(args.out + "/tmp")

#################################################################################### define markers

if ( args.marker == "COI-5P" ):
	marker_list=[ "COI-5P" ]
	coding = True
	if ( args.code == None ):
		print("\n WARNING! a genetic code has to be specified when using a coding marker! \n")
		quit()
if ( args.marker == "COI-3P" ):
	marker_list=[ "COI-3P" ]
	coding = True
	if ( args.code == None ):
		print("\n WARNING! a genetic code has to be specified when using a coding marker! \n")
		quit()
elif ( args.marker == "rbcL" ):
	marker_list=[ "rbcl" , "RBCL" , "Rbcl" , "rbcla" , "rbcl-a" ]
	coding = True
	if ( args.code == None ):
		print("\n WARNING! a genetic code has to be specified when using a coding marker! \n")
		quit()
elif ( args.marker == "matK" ):
	marker_list=[ "matK" , "matk" , "MATK" ]
	coding = True
	if ( args.code == None ):
		print("\n WARNING! a genetic code has to be specified when using a coding marker! \n")
		quit()
elif ( args.marker == "ITS" ):
	marker_list=[ "ITS" , "its" ]
	coding = False
	if ( args.code != None ):
		print("\n WARNING! a genetic code has been specified when using a non-coding marker! \n")
		quit()

if ( args.custom_fas != "" ) and ( args.marker == "coding" ):
	coding = True
	if ( args.code == None ):
		print("\n WARNING! a genetic code has to be specified when using a coding marker! \n")
		quit()
elif ( args.custom_fas != "" ) and ( args.marker == "noncoding" ):
	coding = False
	if ( args.code != None ):
		print("\n WARNING! a genetic code has been specified for a non-coding marker! \n")
		quit()

#################################################################################### big IF

if (args.custom_fas == ""):

#################################################################################### download from BOLD
	
	if ( args.geo == "" ):
		print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t downloading" , args.marker , "sequences for taxa" , args.taxa)
		subprocess.run(["../../bold-cli" , "-taxon" , args.taxa , "-marker" , args.marker , "-output" , "tmp.bold"] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)
		if os.path.getsize("tmp.bold") == 0 :
			print("\n WARNING! The BOLD search for marker" , args.marker , "and taxa" , args.taxa , "returned nothing! \n")
			quit()
	else:
		print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t downloading" , args.marker , "sequences for taxa" , args.taxa, "from" , args.geo)
		subprocess.run(["../../bold-cli" , "-taxon" , args.taxa , "-marker" , args.marker , "-geo" , args.geo, "-output" , "tmp.bold"] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)
		if os.path.getsize("tmp.bold") == 0 :
			print("\n WARNING! The BOLD search for marker" , args.marker , "and taxa" , args.taxa , "in" , args.geo , "returned nothing! \n")
			quit()
						
	num_lines = sum(1 for line in open('tmp.bold', errors='ignore'))
	print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t downloaded" , num_lines , "sequences")

#################################################################################### filter taxonomy and remove gaps
	
	tmp_1_fasta=[]
	
	if (args.filter_tax != ""):
		print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t filtering sequences with" , args.filter_tax , "identification")
	else:
		print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t skipping taxonomic filtering")
	
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

#################################################################################### IF CODING - find longest orfs and length cutoff

if (coding == True):
	
	print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t extracting longest orf - minimum of" , args.filter_len , "amminoacids")
	
	subprocess.run(["getorf" , "-sequence", "tmp2.fna" , "-outseq" , "tmp1.faa" , "-table" , args.code] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)
	
	tmp1_faa=SeqIO.parse("tmp1.faa","fasta")
	tmp2_faa=open("tmp2.faa",'w')
	
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

	tmp3_fna=open("tmp3.fna",'w')

	for record in SeqIO.parse("tmp2.fna", "fasta"):
		if record.id in selected_seqs_header:
			selected_seqs_nt.append(record)
	SeqIO.write(selected_seqs_nt, tmp3_fna, "fasta")   

#################################################################################### IF NON-CODING - length cutoff

if (coding == False):

	tmp3_fna=open("tmp3.fna",'w')

	selected_seqs_nt = list()

	for record in SeqIO.parse("tmp2.fna", "fasta"):
		if len(record.seq) >= int(args.filter_len):
			selected_seqs_nt.append(record)
	SeqIO.write(selected_seqs_nt , tmp3_fna , "fasta")

#if num < 4 :
#	print("\n WARNING! Less than 4 sequences passed the length cutoff! \n" , num)
#	quit()

#################################################################################### IF CODING - align aminoacids and retrotranslate alignment

if (coding == True):

	print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t aligning sequences")
	
	with open('tmp3.faa', 'w') as tmp3_faa, open('tmp2.faa', 'r+') as tmp2_faa:
		os.environ.pop('MAFFT_BINARIES')
		if '*' in tmp2_faa.read():
	 		print("\n WARNING! Stop codons are present in the MSA alignment .. have you selected the right genetic code? \n" , num)
	 		quit()
		
		subprocess.call(["mafft" , "--adjustdirection" , "tmp2.faa"], stdout=tmp3_faa, stderr=subprocess.DEVNULL)

	with open('tmp.aln', 'w') as tmp_aln , open('tmp3.fna', 'r') as tmp3_fna:
		subprocess.run(["pal2nal.pl" , "tmp3.faa" , "tmp3.fna" , "-codontable" , args.code , "-output" , "fasta"], stdout=tmp_aln, stderr=subprocess.DEVNULL)

#################################################################################### IF NON-CODING - align nucleotides
	
if (coding == False):

	print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t aligning sequences")
	with open('tmp.aln', 'w') as tmp_aln, open('tmp3.fna', 'w') as tmp3_fna:
		os.environ.pop('MAFFT_BINARIES')
		subprocess.call(["mafft" , "--adjustdirection" , "tmp3.fna"], stdout=tmp_aln, stderr=subprocess.DEVNULL)

#################################################################################### fragment size selection

print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t selecting optimal marker size")

####################################################################################  IF CODING - write partition file and infer tree or detect custom nwk

if (args.custom_nwk != ""):

	custom_tre = "../../" + args.custom_nwk

	with open(custom_tre, 'r') as tre_file:
		for l in tre_file:
			nwk_sp_list=l.replace(',', ' ').replace('(', '').replace(')', '').replace(';', '').split()

	aln_sp_list = []
	for record in SeqIO.parse("tmp.aln", "fasta"):
		aln_sp_list.append(record.id)

	check =  all(sp in aln_sp_list for sp in nwk_sp_list)

	if  check == False:
		print("\n WARNING! Some species in the phylogeny are not present in the msa! \n")
		quit()
	elif len(nwk_sp_list) == len(aln_sp_list) and check == True:
		print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t custom nwk detected with 1:1 species match to alignmnt - skipping phylogenetic inference")
	elif len(nwk_sp_list) < len(aln_sp_list)  and check == True:
		print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t custom nwk detected with fewer species than alignment - using constrained topology")
	
		if (coding == True):

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
				subprocess.run(["iqtree" , "-s" , "tmp.aln" , "-nt" , str(args.cores) , "-spp" , "tmp.partitions" , "-g" , custom_tre], stdout=tmp_nwk, stderr=subprocess.DEVNULL)
				
			def_tre_file = args.out + ".nwk"
			os.rename('tmp.partitions.treefile' , def_tre_file)
			
			def_aln_file = args.out + ".aln"
			os.rename('tmp.aln' , def_aln_file)
			
		if (coding == False):

			with open('tmp.nwk', 'w') as tmp_nwk :
				subprocess.run(["iqtree" , "-s" , "tmp.aln" , "-nt" , str(args.cores) , "-g" , custom_tre], stdout=tmp_nwk, stderr=subprocess.DEVNULL)
				
			def_tre_file = args.out + ".nwk"
			os.rename('tmp.partitions.treefile' , def_tre_file)
			
			def_aln_file = args.out + ".aln"
			os.rename('tmp.aln' , def_aln_file)

elif (args.custom_nwk == ""):
	
	print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t inferring tree")
	
	if (coding == True):

		tmp_partitions=[]
		filename = "tmp.aln"
		format = "fasta"
		tmp_aln = AlignIO.read("tmp.aln", "fasta")
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

	if (coding == False):

		with open('tmp.nwk', 'w') as tmp_nwk :
			subprocess.run(["iqtree" , "-s" , "tmp.aln" , "-nt" , str(args.cores) ], stdout=tmp_nwk, stderr=subprocess.DEVNULL)
			
		def_tre_file = args.out + ".nwk"
		os.rename('tmp.aln.treefile' , def_tre_file)
		
		def_aln_file = args.out + ".aln"
		os.rename('tmp.aln' , def_aln_file)

################################################################################### ancestral sequences inference

if (coding == True):
	print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t inferring ancestral sequences using genetic code" , args.code)
if (coding == False):
	print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t inferring ancestral sequences")

tmp_ctl=[]

aln_line="seqfile = " + def_aln_file
tmp_ctl.append(aln_line)
tmp_ctl.append("outfile = tmp_baseml.out")
tre_line="treefile = " + def_tre_file
tmp_ctl.append(tre_line)
tmp_ctl.append("noisy = 3")
tmp_ctl.append("verbose = 1")
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
	icode = int(args.code) - 1
	icode_line="icode = " + str(icode)
	tmp_ctl.append(icode_line)
tmp_ctl.append("fix_blength = 2")
tmp_ctl.append("method = 0")

with open('tmp.ctl', 'w') as tmp_ctl_file:
    for line in tmp_ctl:
            tmp_ctl_file.write(line + '\n')
            
with open('tmp.nwk', 'w') as tmp_nwk, open('tmp_baseml.log', 'w') as tmp_baseml_log:
	subprocess.run(["baseml" , "tmp.ctl"], stdout=tmp_baseml_log, stderr=subprocess.DEVNULL)

################################################################################### extract baseml trees

with open("rst", 'r') as rst, open("brlen.phylo.txt", 'w') as brlen_phylo, open("nodes.phylo.txt", 'w') as nodes_phylo:
    brlen = [10]
    nodes = [17]
    for position, line in enumerate(rst):
        if position in brlen:
            brlen_phylo.write(line.replace(" ", ""))
        if position in nodes:
            nodes_phylo.write(line.replace(" ", ""))
            
################################################################################### RRRRRR

print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t finding baits with max distance of", args.distance, "nt/site")

tmp_rscript=[]
tmp_rscript.append("require(geiger)")
tmp_rscript.append("require(ape)")
tmp_rscript.append("require(RRphylo)")
tmp_rscript.append("require(adephylo)")
tmp_rscript.append("require(phytools)")
tmp_rscript.append("library(phangorn)")
tmp_rscript.append("library(data.table)")
tmp_rscript.append("")
tmp_rscript.append("###################################################################### params ####")
tmp_rscript.append("")
param_line="param=" + args.distance
tmp_rscript.append(param_line)
tmp_rscript.append("  grain=param")
tmp_rscript.append("  distance_cutoff=param")
tmp_rscript.append("")
tmp_rscript.append("###################################################################### script 1 ####")
tmp_rscript.append("")
tmp_rscript.append("treebr<-read.newick(\"brlen.phylo.txt\")")
tmp_rscript.append("treend<-read.newick(\"nodes.phylo.txt\")")
tmp_rscript.append("treebr$node.label <- treend$node.label")
tmp_rscript.append("tree <- as.phylo(treebr)")
tmp_rscript.append("")
tmp_rscript.append("up <- sort((distRoot(tree, tips = \"all\", method =\"patristic\")), decreasing = TRUE)[1]")
tmp_rscript.append("up <- up-(up/100)")
tmp_rscript.append("ndf <- NULL")
tmp_rscript.append("ndf_col <- list()")
tmp_rscript.append("")
tmp_rscript.append("root<-length(tree$tip)+1")
tmp_rscript.append("node.height<-matrix(NA,nrow(tree$edge),2)")
tmp_rscript.append("for(i in 1:nrow(tree$edge)){")
tmp_rscript.append("  if(tree$edge[i,1]==root){")
tmp_rscript.append("    node.height[i,1]<-0.0")
tmp_rscript.append("    node.height[i,2]<-tree$edge.length[i]")
tmp_rscript.append("  } else {")
tmp_rscript.append("    node.height[i,1]<-node.height[match(tree$edge[i,1],tree$edge[,2]),2]")
tmp_rscript.append("    node.height[i,2]<-node.height[i,1]+tree$edge.length[i]")
tmp_rscript.append("  }")
tmp_rscript.append("}")
tmp_rscript.append("")
tmp_rscript.append("###################################################################### script 2 ####")
tmp_rscript.append("")
tmp_rscript.append("for (slice in seq(up, 0, by = -grain)){")
tmp_rscript.append("  edgs<-which(node.height[,2]>slice&node.height[,1]<slice)")
tmp_rscript.append("  nods<-tree$edge[edgs,2]; nods<-nods[nods>length(tree$tip)]")
tmp_rscript.append("  for (n in 1:length(nods)) {                                                                 # for each of those nodes")
tmp_rscript.append("    nodetips <- NULL                                                                          # empty nodetips")
tmp_rscript.append("    if (! is.na(nods[n])) {")
tmp_rscript.append("      distances <- distNodes(tree,node=nods[n], clus=4)                                       # calculate the distances with all edges")
tmp_rscript.append("      for(i in 1:nrow(distances)) {                                                           # for each row of the distance df")
tmp_rscript.append("        t <- row.names(distances)[i]                                                          #  ")
tmp_rscript.append("        if (is.element(t,tree$tip.label) == TRUE) {                                           # if the distance is with a terminal")
tmp_rscript.append("          if (distances[i + (nrow(distances))] < distance_cutoff) {                           # and its distance is less than a cutoff (0.5)")
tmp_rscript.append("            nodetips <- append(nodetips, t) ")
tmp_rscript.append("          }")
tmp_rscript.append("        }")
tmp_rscript.append("      }")
tmp_rscript.append("      print(paste0(\"nodes \",nods[n],\" at cutoff \",round(slice, digits=3),\" - includes \", length(nodetips), \" species\"))")
tmp_rscript.append("      print(nodetips)")
tmp_rscript.append("      if (! is.null(nodetips)) {")
tmp_rscript.append("        for (d in 1:length(nodetips)) { ndf[ndf == nodetips[d]] <- NA }                       # turn to NA all tips previously")
tmp_rscript.append("        length(nodetips) <- max(length(tree$tip.label))")
tmp_rscript.append("        #    names(nodetips) <- nods[n]")
tmp_rscript.append("        ndf <- cbind(ndf, lapply(nodetips, unlist))")
tmp_rscript.append("        ndf_col <- append(ndf_col, nods[n])")
tmp_rscript.append("      }")
tmp_rscript.append("    }")
tmp_rscript.append("  }")
tmp_rscript.append("}")
tmp_rscript.append("")
tmp_rscript.append("###################################################################### polish ####")
tmp_rscript.append("")
tmp_rscript.append("ndf <- as.data.frame(ndf)")
tmp_rscript.append("colnames(ndf) <- ndf_col")
tmp_rscript.append("ndf <- ndf[,colSums(is.na(ndf))<nrow(ndf)]")
tmp_rscript.append("")
tmp_rscript.append("for (m in 1:length(colnames(ndf))) { t=as.integer(colnames(ndf)[m])-as.integer(tree$Nnode); if (t > 0) { print((tree$node.label[t])) } }")
tmp_rscript.append("")
tmp_rscript.append("for (m in 1:length(colnames(ndf))) {print(as.integer(colnames(ndf)[m]))}")
tmp_rscript.append("")
tmp_rscript.append("plot.phylo(tree,cex=1, type=\"unrooted\", use.edge.length=TRUE); nodelabels()")
tmp_rscript.append("for (y in 1:length(colnames(ndf))) {nodelabels("", as.integer(colnames(ndf)[y]), frame = \"c\", bg = \"yellow\")}")
tmp_rscript.append("")
tmp_rscript.append("sum(colSums(!is.na(ndf)))")
tmp_rscript.append("fwrite(list(colnames(ndf)), file = \"nodes.csv\")")

with open('tmp.tmp_rscript', 'w') as tmp_rscript_file:
    for line in tmp_rscript:
            tmp_rscript_file.write(line + '\n')

with open('tmp_rscript_log', 'w') as tmp_rscript_log:
	subprocess.run(["Rscript" , "tmp.tmp_rscript"], stdout=tmp_rscript_log, stderr=subprocess.DEVNULL)

################################################################################### extract probes

with open("nodes.csv", 'r') as nodes:
	for linenodes in nodes:
		pattern=("node #"+linenodes).strip()
		with open("rst", 'r') as rst, open ("baits", 'a') as baits:
			for line in rst:
				if re.search(pattern, line):
					baits.write((line.replace(" ", "")).replace("node#[0-9]*", ">bait \n"))

################################################################################### clean

os.chdir('..')
shutil.copy('tmp/baits', '.')

if args.verbose == False :
	pass
else :
	files = glob.glob("tmp/*")
	for f in files:
		os.remove(f)
	os.rmdir("tmp")

print("\n" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), ":\t finish \n")