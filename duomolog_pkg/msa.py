from Bio.Align.Applications import MafftCommandline
from io import StringIO
from Bio import AlignIO

def run_mafft(inFile):
	mafft_cline = MafftCommandline(input = inFile)
	stdout, stderr = mafft_cline()
	alignment = AlignIO.read(StringIO(stdout), "fasta")
	
	# print(align)
	return alignment
	