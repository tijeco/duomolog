# AUTOGENERATED! DO NOT EDIT! File to edit: notebooks/01_blast.ipynb (unless otherwise specified).

__all__ = ['run_blast', 'load_blast']

# Cell
from io import StringIO
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
import pandas as pd

# Cell
def run_blast(dbFile, queryFile):
	makeblastdbCMD = NcbimakeblastdbCommandline(dbtype="prot",
    	input_file=dbFile)
	blastpCMD = NcbiblastpCommandline(query=queryFile,
		db=dbFile, evalue=0.001, max_target_seqs=1, outfmt=6)
	print(makeblastdbCMD)
	print(blastpCMD)
	makeblastdbOUT, makeblastdbERR = makeblastdbCMD()
	
	blastpOUT, blastpERR = blastpCMD()
	blastpDF = pd.read_csv(StringIO(blastpOUT),sep="\t", names = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
	blastpDF["sseqid"].astype(str)
	blastpDF["qseqid"].astype(str)
	
	return blastpDF
		

# Cell
def load_blast(blastout,inHeaders,hasHeaders=False):
	if hasHeaders:
		blastpDF = pd.read_csv(blastout, sep = "\t")
	else:
		blastpDF = pd.read_csv(blastout, sep = "\t", names = ["qseqid","sseqid","pident","length","mismatch",			"gapopen","qstart","qend","sstart","send","evalue","bitscore"])
	blastpDF["sseqid"].astype(str)
	blastpDF["qseqid"].astype(str)

	blastHeaders = set(blastpDF["sseqid"])
	blastHeaders_union_inHeaders = blastHeaders | set(inHeaders)
	print(blastHeaders)
	print(set(inHeaders))
	if blastHeaders_union_inHeaders  != set(inHeaders):
		sys.exit("the provided blast file contains headers not found in the input file")
	
	return blastpDF