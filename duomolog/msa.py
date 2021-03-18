# AUTOGENERATED! DO NOT EDIT! File to edit: notebooks/02_msa.ipynb (unless otherwise specified).

__all__ = ['run_mafft', 'writeAln', 'fa2sto']

# Cell
from Bio.Align.Applications import MafftCommandline
from io import StringIO
from Bio import AlignIO

# Cell
def run_mafft(fastaFile,format="fasta"):
	mafft_cline = MafftCommandline(input = fastaFile)
	print(mafft_cline)
	stdout, stderr = mafft_cline()
	alignment = AlignIO.read(StringIO(stdout), format)
	alignmentFile = writeAln(alignment, fastaFile,format="fasta")
		
	# print(align)
	return alignmentFile, alignment

# Cell
def writeAln(aln, fastaFile,format="fasta"):
    alnFile = fastaFile + ".aln"
    AlignIO.write(aln,alnFile,format)
    return alnFile

# Cell
def fa2sto(faFile):
    fa = AlignIO.read(faFile,"fasta")
    stoFile = faFile + ".sto"
    AlignIO.write(fa,stoFile,"stockholm")


    return stoFile