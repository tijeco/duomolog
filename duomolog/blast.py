from io import StringIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
import pandas as pd 

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
	
	return blastpDF

def load_blast(blastout,inHeaders):
	blastpDF = pd.read_csv(blastout, sep = "\t", names = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
	blastHeaders = set(blastpDF["sseqid"].astype(str))
	blastHeaders_union_inHeaders = blastHeaders | set(inHeaders)
	if blastHeaders_union_inHeaders  != set(inHeaders):
		sys.exit("the provided blast file contains headers not found in the input file")
	
	return blastpDF
		