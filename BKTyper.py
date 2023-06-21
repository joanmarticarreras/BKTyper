#!/usr/bin/env python3

# Import packages
from Bio import SeqIO
from Bio import AlignIO
from Bio.Emboss.Applications import NeedleCommandline
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from io import StringIO
import Bio
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors
from Bio.Align.Applications import MafftCommandline
import os
import sys
import re
import regex
import numpy as np
import pandas as pn
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import pylab
import subprocess
from Bio.Application import _Option, _Switch, AbstractCommandline
from Bio.Phylo.Applications import PhymlCommandline
from Bio.SeqRecord import SeqRecord
from PyPDF2 import PdfFileMerger
import shutil

# Set globals
pn.set_option('mode.chained_assignment', None)


###################
#######USAGE#######
###################

# $ python3 BKTyper.py <input-sequences> <mode: VP1, NCCR, complete>

###################

# Functions

def NCCR_classification(blast_out):
	# NCCR Block classification with BLAST
	array_block = blast_out.split("\n")
	array_block = (pn.DataFrame(array_block[:-1]))
	dataframe = array_block[0].str.split("\t",expand=True)
	dataframe.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','score']
	dataframe.qstart = dataframe.qstart.astype(int)
	dataframe.qend = dataframe.qend.astype(int)
	dataframe_sorted = dataframe.sort_values('qstart')
	dataframe_sorted_nr = dataframe_sorted.drop_duplicates(subset=["qstart","qend"], keep='last')
	dataframe_sorted_nr.qstart = dataframe_sorted_nr.qstart.astype(str)
	dataframe_sorted_nr.qend = dataframe_sorted_nr.qend.astype(str)
	dataframe_sorted_nr['NCCR_positions'] = dataframe_sorted_nr.qstart + '-' + dataframe_sorted_nr.sseqid + '-' + dataframe_sorted_nr.qend
	NCCR = ''.join(dataframe_sorted_nr['sseqid'].tolist())
	NCCR_complex = '|'.join(dataframe_sorted_nr['NCCR_positions'].tolist())
	return (NCCR,NCCR_complex)


def motif_finder(sequence,motif_file,NCCR):
	motifs = []
	motif_file = motif_file.read().splitlines()
	NCCR_length = int(NCCR[-3:])
	for i in motif_file:
		(id,nucleotides,length) = i.split("\t")
		motif = (re.finditer(nucleotides,sequence[1:NCCR_length]))
		for m in motif:
			start = m.start()+1
			end = m.end()+1
			motif = (id,start,end)
			motifs.append(motif)
	blocks = NCCR.split("|")
	for i in blocks:
		(Start,Block,End) = tuple(i.split("-"))
		block = (Block,int(Start),int(End))
		motifs.append(block)
	df = pn.DataFrame(motifs,columns = ["Motif","Start","End"])
	df_ordered = df.sort_values(by=["Start"])
	return(df_ordered)


def write_visuals(seq_name,df,seq_length,results_pdf):
	gd_diagram = GenomeDiagram.Diagram(seq_name, track_size=1)
	new_row = pn.DataFrame({"Motif":seq_name,"Start":1,"End":1},index=[0])
	df = pn.concat([new_row, df]).reset_index(drop = True)
	NCCR = GenomeDiagram.FeatureSet()
	for index, row in df.iterrows():
		(motif_name,start_motif,end_motif) = (row["Motif"],row["Start"],row["End"])
		cols = [motif_name,start_motif,end_motif]
		if index == 0:
			block = SeqFeature(FeatureLocation(int(cols[1]), int(cols[2]), strand=-1), type="blocks", id=motif_name)
			NCCR.add_feature(block,color=colors.HexColor("#8DD35F"), name=motif_name, label=True, label_size=8, label_position="middle", label_angle=180)
		else:
			if motif_name.islower():
				motif = SeqFeature(FeatureLocation(int(cols[1]), int(cols[2]), strand=+1), type="motifs",id=motif_name)
				NCCR.add_feature(motif, color=colors.HexColor("#8DD35F"), name=motif_name, label=True, label_size=10, label_position="left", label_angle=90)
			else:
				block = SeqFeature(FeatureLocation(int(cols[1]), int(cols[2]), strand=-1), type="blocks", id=motif_name)
				if motif_name == "O":
					NCCR.add_feature(block,color=colors.HexColor("#ffc69e"), name=motif_name, label=True, label_size=10, label_position="center", label_angle=180)
				elif motif_name == "P":
					NCCR.add_feature(block,color=colors.HexColor("#fff6d4"), name=motif_name, label=True, label_size=10, label_position="middle", label_angle=180)
				elif motif_name == "Q":
					NCCR.add_feature(block,color=colors.HexColor("#f6f9eb"), name=motif_name, label=True, label_size=10, label_position="middle", label_angle=180)
				elif motif_name == "R":
					NCCR.add_feature(block,color=colors.HexColor("#ebf9f6"), name=motif_name, label=True, label_size=10, label_position="middle", label_angle=180)
				elif motif_name == "S":
					NCCR.add_feature(block,color=colors.HexColor("#f9ebf6"), name=motif_name, label=True, label_size=10, label_position="middle", label_angle=180)
				else:
					NCCR.add_feature(block,color=colors.HexColor("#C8C4B7"), name=motif_name, label=True, label_size=20, label_position="right", label_angle=180)
	NCCR_track =  GenomeDiagram.Track(name="Annotated Features", height=0.3)
	NCCR_track.add_set(NCCR)
	gd_diagram.add_track(NCCR_track,3)
	seq_length = int(cols[2])
	rows = max(2, int(round(seq_length / 100)))
	gd_diagram.draw(format='linear', tracklines=0, pagesize='A4', orientation = 'landscape', fragments=4, start=1, end=int(seq_length))

	pdf_filepath = os.path.join('results','{}.pdf'.format(seq_name))
	gd_diagram.write(pdf_filepath, 'PDF', dpi=300)
	results_pdf.append(pdf_filepath)

def VP1_classification(aln):
	# Decision tree for Polyoma BK subtyping based on VP1 polymorphisms (Morel et al Journal of Clinical Microbiology 2017, https://doi.org/10.11128/JCM.01180-16)
	subgroup = "NA"
	coordinates = {426:'aaaacctat',513:'aaagtac',456:'ctttgctg',465:'aggtggagaa',444:'taatttccacttctttg',495:'gctaatgaattacag',433:'attcaaggcagtaattt',414:'gcatggtggaggaaa'}
	for position in coordinates:
		for i in re.finditer(coordinates[position],str(aln[1,:].seq)):
			coordinates[position] = i.start()

	if aln[:,coordinates[426]][0].upper() == "G".upper():
    		subgroup = "Ic"
	elif aln[:,coordinates[513]][0].upper() == "C".upper():
    		subgroup = "Ib-1"
	elif aln[:,coordinates[426]][0].upper() == "A".upper():
   		subgroup = "Ia"
	elif aln[:,coordinates[456]][0].upper() == "T".upper() and aln[:,coordinates[465]][0].upper() == "T".upper():
   		subgroup = "Ib-2"
	elif aln[:,coordinates[426]][0].upper() == "T".upper():
   		subgroup = "III"
	elif aln[:,coordinates[444]][0].upper() == "T".upper():
   		subgroup = "II"
	elif aln[:,coordinates[495]][0].upper() == "G".upper():
   		subgroup = "IVa-1"
	elif aln[:,coordinates[495]][0].upper() == "C".upper():
   		subgroup = "IVa-2"
	elif aln[:,coordinates[433]][0].upper() == "G".upper():
   		subgroup = "IVc-1"
	elif aln[:,coordinates[414]][0].upper() == "A".upper():
   		subgroup = "IVc-2"
	else:
   		subgroup = "IVb-1,2"

	subgroup_detail = aln[:,coordinates[426]][0].upper()+aln[:,coordinates[513]][0].upper()+aln[:,coordinates[456]][0].upper()+aln[:,coordinates[465]][0].upper()+aln[:,coordinates[444]][0].upper()+aln[:,coordinates[495]][0].upper()+aln[:,coordinates[433]][0].upper()+aln[:,coordinates[414]][0].upper()
	return(subgroup,subgroup_detail)

def VP1_tree(seq_name,seq,db):
	sequences = ""
	sequences = sequences.join([db,"\n",">",seq_name,"\n",seq,"\n"])
	sequences_file = open("query_mafft_vp1.fasta", "w")
	sequences_file.write(sequences)
	sequences_file.close()
	mafft_cline = MafftCommandline(input='query_mafft_vp1.fasta')
	(mafft_stdout,mafft_stderr) = mafft_cline()
	with open("aligned.fasta", "w") as handle:
		handle.write(mafft_stdout)
	mafft_align = AlignIO.read(StringIO(mafft_stdout), "fasta")
	VP1_alignment = trim_alignment(mafft_align)
	BKTGR_alignment = VP1_alignment[:,425:512]
	AlignIO.write(BKTGR_alignment, "BKTGR_alignment.aln", "clustal")
	os.system("iqtree -s BKTGR_alignment.aln -bb 10000 -redo -m TEST -pre model_testing")
	tree = Phylo.read("model_testing.contree", "newick")
	for node in tree.get_nonterminals():
		node.name = None
	tree.ladderize()
	tree.root_at_midpoint()
	matplotlib.rc('font', size=5)
	fig = pylab.figure(figsize=(10, 30))
	axes = fig.add_subplot(1, 1, 1)
	morel_vp1=[]
	with open("source/VP1_BKTyper_MLtree_list.txt") as vp1_typing:
		next(vp1_typing)
		for i in vp1_typing:
			(q,vp1) = tuple(i.split("\t"))
			line = ("num",q,"s","nccr","dnccr",vp1,"dvp1")
			morel_vp1.append(line)
	morel_vp1.append(("num",seq_name,"s","nccr","dnccr","Query","dvp1"))
	df = pn.DataFrame(morel_vp1,columns=["Num","Query","Strand","NCCR","Detail_NCCR","VP1_Subtype","VP1_Subtype_Detailed"])
	df = df.replace('\n','', regex=True)
	label_color = {}
	color_legend = {"Query":"#FF0000","Ia":"#FF8080","Ib-1":"#FF9955","Ib-2":"#AC9D93","Ic":"#FFE680","II":"#37C871","III":"#8DD35F","IVa-1":"#93AC9D","IVa-2":"#80E5FF","IVb-1":"#8787DE","IVb-2":"#8787DE","IVc-1":"#E580FF","IVc-2":"#DE87AA"}
	for index, row in df.iterrows():
		(name,type) = row["Query"],row["VP1_Subtype"]
		color = color_legend[type]
		label_color[name] = color
	pyplot.rcParams.update({'font.size': 18})
	Phylo.draw(tree, do_show=False,axes=axes,label_colors=label_color,show_confidence=False)
	pylab.title(seq_name,loc="center")
	pylab.axis("off")
	pylab.savefig("iqtree_GTR_fast.svg",format='svg', papertype = "a5", transparent = True, dpi=300)

def trim_alignment(alignment):
	start_VP1 = 0
	end_VP1 = 1
	for col in range(alignment.get_alignment_length()):
		if not "-" in alignment[:, col]:
			start_VP1 = col
			break
		else:
			pass
	for col in reversed(range(alignment.get_alignment_length())):
		if not "-" in alignment[:, col]:
			end_VP1 = col
			break
		else:
			pass
	trim_alignment = alignment[:, start_VP1:end_VP1]
	AlignIO.write(trim_alignment, "test.aln", "clustal")
	return(trim_alignment)

# Create outputs
end_table = pn.DataFrame(columns=['Query','Strand','NCCR','Detail_NCCR','VP1_Subtype','VP1_Subtype_Detailed'])
convert_fasta = []
results_pdf = []
os.mkdir("results")
# Open (multi)fasta
for sequence in SeqIO.parse(sys.argv[1],"fasta"):
	sequence.seq = sequence.seq.upper()
	SeqIO.write(sequence,"target_sequence","fasta")
# Detect and correct strandess of the sequence with BLAST
	strand = NcbiblastnCommandline(query="target_sequence",subject="source/dunlop.fasta", outfmt=5, max_hsps=1)()[0]
	blast_result_record = NCBIXML.read(StringIO(strand))
	sign = ""
	for description in blast_result_record.descriptions:
		if (description.e < 0.000001):
			pass
		else:
			print(sequence.id+" Input sequence is not similar enough to polyoma BK")
			continue

	for alignment in blast_result_record.alignments:
		for hsp in alignment.hsps:
			if ( (hsp.sbjct_end - hsp.sbjct_start) > 0 ):
				sign = "+"
			else:
				sequence = sequence.reverse_complement()
				sign = "-"
				SeqIO.write(sequence,"target_sequence","fasta")


# Alignments
	if sys.argv[2] == 'complete':
		if re.search('TTTTGC(.AAAA|A.AA|AA.A|AAA.)',str(sequence.seq)) is not None:
			if re.search('TTTTGC(.AAAA|A.AA|AA.A|AAA.)',str(sequence.seq)).start() > 0:
				seq_pieces = re.split('(TTTTTGC[.AAAA|A.AA|AA.A|AAA.])',str(sequence.seq))
				seq_pieces_fix = seq_pieces[1]
				sequence.seq = Bio.Seq.Seq(seq_pieces_fix[1:] + seq_pieces[2] + seq_pieces[0] + "T")
				SeqIO.write(sequence,"target_sequence","fasta") # Output the input sequence restructured as Dunlop reference
		else:
			print(sequence.id+" lacks the origin of replication")
			continue

	# NCCR BLAST
	block =  NcbiblastnCommandline(query="target_sequence", subject="source/NCCR_BKTyper.fasta", outfmt=6, word_size=12, perc_identity=75, evalue=0.05)()[0] ###
	# VP1 Needleman and Wunch
	a = NeedleCommandline(asequence="target_sequence", \
				bsequence="source/VP1_Dunlop.fasta", \
				gapopen=10, \
				gapextend=0.5, \
				outfile="needle_fname")

	a() # execute the alignment
	# Export the alignment back to Python
	VP1_alignment = AlignIO.read("needle_fname", "emboss")

# Call functions based on mode
	NCCR=NCCR_complex=subgroup=subgroup_detail = 'NA' # definition of table objects
	motif_list=(open("source/motif_list.txt","r"))
	vp1_db_file=(open("source/VP1_BKTyper_MLtree_list.fasta","r"))
	vp1_db = vp1_db_file.read()

	if sys.argv[2] == 'VP1':
		(subgroup,subgroup_detail) = VP1_classification(VP1_alignment)
		if sys.argv[3] == 'Tree':
			VP1_tree(str(sequence.id),str(sequence.seq),vp1_db)
		else:
			pass
	elif sys.argv[2] == 'NCCR':
		(NCCR,NCCR_complex) = NCCR_classification(block)
		regions = motif_finder(str(sequence.seq),motif_list,NCCR_complex)
		length=len(sequence.seq)
		write_visuals(sequence.id,regions,length,results_pdf)
	elif sys.argv[2] == 'complete':
		block =  NcbiblastnCommandline(query="target_sequence", subject="source/NCCR_BKTyper.fasta", outfmt=6, word_size=12, perc_identity=75, evalue=0.05)()[0] ### # if complete genome is provided, NCCR numbering should be based on the re-ordered sequence
		(subgroup,subgroup_detail) = VP1_classification(VP1_alignment)
		if sys.argv[3] == 'Tree':
			VP1_tree(str(sequence.id),str(sequence.seq),vp1_db)
		else:
			pass
		(NCCR,NCCR_complex) = NCCR_classification(block)
		regions = motif_finder(str(sequence.seq),motif_list,NCCR_complex)
		length=len(sequence.seq)
		write_visuals(sequence.id,regions,length,results_pdf)
	else:
		sys.exit("Select analysis mode: NCCR, VP1 or complete" + "\n")
	# Append results to the results table and (multi)fasta file
	end_table = end_table.append({'Query':sequence.id,'Strand':sign,'NCCR':NCCR,'Detail_NCCR':NCCR_complex,'VP1_Subtype':subgroup,'VP1_Subtype_Detailed':subgroup_detail},ignore_index=True)
	convert_fasta.append(sequence)

end_table.to_csv(r'results.csv',sep="\t")
SeqIO.write(convert_fasta,"results.fasta","fasta")

merger = PdfFileMerger()
for pdf in results_pdf:
	merger.append(pdf)
merger.write("results.pdf")
merger.close()

# Delete intermediate files
os.remove("target_sequence")
os.remove("needle_fname")
shutil.rmtree("results/")
