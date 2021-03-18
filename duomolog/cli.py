# AUTOGENERATED! DO NOT EDIT! File to edit: notebooks/00_cli.ipynb (unless otherwise specified).

__all__ = ['Display', 'main', 'CLI', 'Duo', 'writeOut', 'nrPepDB', 'pepDB2Fa', 'subPepDB', 'intersectInQuery']

# Cell
import argparse, os, sys
import fire
from fire import core
from Bio import SeqIO
import pandas as pd
from pathlib import Path

# Cell
from duomolog import blast
from duomolog import msa
from duomolog import hmmer

# Cell
def Display(lines, out):
    text = "\n".join(lines) + "\n"
    out.write(text)
def main():
    commands = CLI()
    core.Display = Display
    fire.Fire(commands, name="duomolog")
    # print(say_hello("Jeremy")=="Hello Jeremy!")

# Cell
class CLI(object):
    def __init__(self):
        self.blastOut = None
        self.hmmerOut = None
    def blast_v_hmmer(self,
		inFile,
		queryFile,
		summaryOut="duomolog_out/summary_out.txt",
		blastFile=None,
        alnFile=None,
		hmmFile=None,
        intersectOnly=False):
        os.makedirs(os.path.dirname(summaryOut), exist_ok=True)
        Path(summaryOut).touch() # should at worst guarantee output is generated
        outDir = '/'.join(summaryOut.split("/")[:-1]) +"/"
        in_nrPepDB = nrPepDB(inFile)
        in_nrPepDB.to_csv(inFile+".nr.csv")
        in_nrFile = outDir + "in_nr.fasta"
        pepDB2Fa(in_nrPepDB,in_nrFile)

        query_nrPepDB = nrPepDB(queryFile)
        query_nrPepDB.to_csv(queryFile+".nr.csv")
        query_novelDB, query_shared_with_inputDB = intersectInQuery(in_nrPepDB,query_nrPepDB)
        query_novel_file = outDir + "query_novel.fasta"
        query_shared_with_input_File = outDir + "query_shared_with_input.fasta"

        if query_shared_with_inputDB != None:
            pepDB2Fa(query_shared_with_inputDB,query_shared_with_input_File)
        pepDB2Fa(query_novelDB,query_novel_file)

        if blastFile == None:
            self.blastOut = blast.run_blast(in_nrFile, query_novel_file)
        else:
            inHeaders = in_nrPepDB["1stHeader"]
            self.blastOut = blast.load_blast(blastFile,inHeaders)

        if hmmFile == None:
            alnFile, alnOut = msa.run_mafft(in_nrFile)
        self.hmmerOut = hmmer.run_hmmsearch(query_novel_file, alnFile, hmmFile)

        blast_headers = set(self.blastOut["qseqid"])
        hmmer_headers = set([hit.name.decode() for hit in self.hmmerOut])
        blast_hmmer_duo = Duo("blast",blast_headers,"hmmer",hmmer_headers)
        writeOut(outDir,summaryOut,blast_hmmer_duo,query_novelDB,intersectOnly)


# Cell
class Duo:
	def __init__(self, leftName, left, rightName, right):
		self.left_union_rightName = leftName + "_union_" + rightName
		self.left_not_rightName = leftName + "_not_" + rightName
		self.left_intersect_rightName = leftName + "_intersect_" + rightName
		self.right_not_leftName = rightName + "_not_" + leftName

		self.subsets = {
			leftName: left,
			rightName: right,
			self.left_union_rightName: left | right,
			self.left_not_rightName: left - right,
			self.left_intersect_rightName: left & right,
			self.right_not_leftName: right - left
		}

	def dropEmpty(self):
		emptySubsets = []
		for subset in self.subsets:
			if self.subsets[subset] == set():
				print(subset, "is empty, removing it")
				emptySubsets.append(subset)
		for emptySubset in emptySubsets:
			del self.subsets[emptySubset]

	def dropRedundant(self):
		nrSubSets = []
		redundantSubsets = []
		for subset in self.subsets:
			if self.subsets[subset] not in nrSubSets:
				nrSubSets.append(self.subsets[subset])
			else:
				print(subset, "is redundant, removing it")
				redundantSubsets.append(subset)
		for redundantSubset in redundantSubsets:
			del self.subsets[redundantSubset]
	def intersectOnly(self):
		self.subsets = {self.left_intersect_rightName:self.subsets[self.left_intersect_rightName]}


# Cell
def writeOut(outDir,summaryOut,duo,queryDB,intersectOnly):
    if intersectOnly:
        duo.intersectOnly()
    else:
        duo.dropRedundant()
    duo.dropEmpty()
    if bool(duo):
        with open(summaryOut, "w") as out:
            for subset in duo.subsets:
                subsetDB = subPepDB(queryDB, duo.subsets[subset])
                pepDB2Fa(subsetDB, outDir + f"{subset}.fa")
                for header in duo.subsets[subset]:
                    out.write(f"{header}\t{subset}")


# Cell
def nrPepDB(faFile):
    seqDict = {}
    outDict = {"1stHeader":[],  "allHeaders":[],"record":[],"seq":[]}
    with open(faFile) as fa:
        for record in SeqIO.parse(fa, "fasta"):
            if record.seq not in seqDict:
                seqDict[record.seq] = [record]
            else:
                print(f"WARNING: {record.id} contains a redundant sequence")
                seqDict[record.seq].append(record)

    for seq in seqDict:
        h1 = seqDict[seq][0].id
        # allHs = ';'.join(seqDict[seq])
        allHs = ';'.join([r.id for r in seqDict[seq]])
        outDict["1stHeader"].append(h1)
        outDict["allHeaders"].append(allHs)
        outDict["record"].append(seqDict[seq][0])
        outDict["seq"].append(str(seq))


    outDB = pd.DataFrame.from_dict(outDict)
    outDB["1stHeader"].astype(str)
    return outDB





# Cell
def pepDB2Fa(pepDB, faFile,format="fasta"):
    SeqIO.write(pepDB["record"], faFile,format)


# Cell
def subPepDB(pepDB, headers):
    return pepDB[pepDB["1stHeader"].isin(headers)]


# Cell
def intersectInQuery(inDB,queryDB):
    inSeqs = set(inDB["seq"])
    querySeqs = set(queryDB["seq"])
    # print(inSeqs)
    # print(querySeqs)
    sharedSeqs = inSeqs & querySeqs

    if sharedSeqs == set():
        # titanic[titanic["Age"] > 35]
        outDB = queryDB
        sharedDB = None
    else:
        # print(sharedSeqs)
        sharedDB = queryDB[queryDB["seq"].isin(sharedSeqs)]
        outDB = queryDB[~queryDB["seq"].isin(sharedSeqs)]

    return outDB, sharedDB