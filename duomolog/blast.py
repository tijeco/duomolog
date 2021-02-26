from io import StringIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
import pandas as pd

def run_blast(dbFile, queryFile):
	makeblastdbCMD = NcbimakeblastdbCommandline(dbtype="prot",
    	input_file=dbFile)
	makeblastdbOUT, makeblastdbERR = makeblastdbCMD()
	blastpCMD = NcbiblastpCommandline(query=queryFile, 
		db=dbFile, evalue=0.001, max_target_seqs=1, outfmt=6)
	blastpOUT, blastpERR = blastpCMD()
	blasdpDF = pd.read_csv(StringIO(blastpOUT),sep="\t"
			names = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
	
	return blasdpDF



 
 
 
 
 
 
 
 
 
 


