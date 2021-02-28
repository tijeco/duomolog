import sys
import argparse, os
from Bio import AlignIO
from Bio import SeqIO

from duomolog import blast
from duomolog import msa
from duomolog import hmmer
from duomolog import resultSubset




class ParseCommands(object):
	
	def __init__(self):

		parser = argparse.ArgumentParser(
			description="duomolog",
			usage='''duomolog <command> [<args>] <queryPepFasta>
			The most commonly used duomolog commands are:
				blast_v_hmmer
				hmm_v_hmm''')
		parser.add_argument("command", help="Subcommand to run")
		parser.add_argument('-format',
                    choices=["clustal",
							"emboss",
							"fasta",
							"fasta-m10",
							"ig",
							"maf",
							"mauve",
							"msf",
							"nexus",
							"phylip",
							"phylip-sequential",
							"phylip-relaxed",
							"stockholm"],
					default="fasta",
                    help='Special testing value')
		
		parser.add_argument("-outdir","-o",default="duomolog_out")
        
		# parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
		args = parser.parse_args(sys.argv[1:2])
		self.base_args = parser.parse_args(sys.argv[1:2])
		
		if not hasattr(self, args.command):
			print("Unrecognized command")
			parser.print_help()
			exit(1)
		self.args = parser.parse_args(sys.argv[1:2])
        # use dispatch pattern to invoke method with same name
		getattr(self, args.command)()

	def blast_v_hmmer(self):
		parser = argparse.ArgumentParser(
			description="Runs blast and hmmer")
		
		parser.add_argument("--input","-i", type=argparse.FileType('r'),required=True)
		parser.add_argument("--query","-q",type=argparse.FileType('r'), 
			help="FASTA formatted file containing database of peptides to be searched")
		args = parser.parse_args(sys.argv[2:])

		self.args = args
		return(self.args)
	def hmm_v_hmm(self):
		parser = argparse.ArgumentParser(
			description="Runs hmmer on two separate alignments")
		parser.add_argument("-aln1","-a1", type=argparse.FileType('r'),required=True)
		parser.add_argument("-aln2","-a2", type=argparse.FileType('r'),required=True)
		parser.add_argument("--query","-q",type=argparse.FileType('r'), 
			help="FASTA formatted file containing database of peptides to be searched")
		args = parser.parse_args(sys.argv[2:])

		self.args = args
		return(self.args)
		


def notRedundant(sequences):
	seqs = {}
	for id, record in sequences.items():
		seq = record.seq
		if seq not in seqs:
			seqs[seq] = id
		else:
			seqs = {}
			break
	return seqs

def pairRedundant(input,query):
	nr_input = notRedundant(input)
	nr_query = notRedundant(query)
	
	if nr_input == False:
		sys.exit("input alignment contains redundant sequences, please fix and try again")
	if nr_query == False:
		sys.exit("input query contains redundant sequences, please fix and try again")
	
	input_seqs = set(nr_input.keys())
	query_seqs = set(nr_query.keys())
	shared_seqs = input_seqs & query_seqs
	shared_seqs_query_ids = [nr_query[i] for i in shared_seqs]
	# print(shared_seqs)
	# print(len(input_seqs),len(query_seqs),len(shared_seqs))
	return shared_seqs_query_ids
def checkDir(dirName):
	if not os.path.exists(dirName):
		os.makedirs(dirName)
	elif len(os.listdir(dirName)) != 0:
		print("WARNING: provided output directory is not empty")


def main():
	cli = ParseCommands()

	duomolog_args = cli.args
	subcommand = cli.base_args.command 
	alignment_format = cli.base_args.format
	outdir = cli.base_args.outdir
	checkDir(outdir)

	query_file = duomolog_args.query.name
	querySeq = SeqIO.index(query_file, "fasta")
	# print(querySeq.items())

	if subcommand == "blast_v_hmmer":
		input_file = duomolog_args.input.name
		# alignment = AlignIO.read(open(alignment_file), alignment_format)
		inputSeq = SeqIO.index(input_file, "fasta")
		sharedIDs = pairRedundant(inputSeq,querySeq)
		if sharedIDs:
			print("Error: the following sequences were found in the query that are already in the input alignment, please remove and try again")
			for h in sharedIDs:
				print(h)
			sys.exit()

		blast_results = blast.run_blast(input_file,query_file)
		
		# print(blast_results)
		alignment = msa.run_mafft(input_file)
		alignment_file = input_file + ".mafft.aln"
		alignmentOut = SeqIO.write(alignment, alignment_file, alignment_format)
		hmmer_results = hmmer.run_hmmer(alignment_file, query_file)

		blast_headers = set(blast_results["qseqid"].astype(str))
		hmmer_headers = set([hit.name.decode() for hit in hmmer_results])
		blast_hmmer_subset = resultSubset.Duo("blast",blast_headers,"hmmer",hmmer_headers)
		blast_hmmer_subset.dropEmpty()
		blast_hmmer_subset.dropRedundant()

		# print(querySeq)
		# print(querySeq["119355"])
		blast_hmmer_subsetSeqs = resultSubset.seqSubSet(querySeq,blast_hmmer_subset.subsets)
		with open(outdir +"/duomolog_results.txt", "w") as summary_out:
			for subset in blast_hmmer_subsetSeqs:
				with open(outdir +"/"+  subset + ".fa","w") as seq_out:
					for input_header in inputSeq:
						record = inputSeq[input_header]
						seq_out.write(record.format("fasta"))
					for header in blast_hmmer_subsetSeqs[subset]:
						record = blast_hmmer_subsetSeqs[subset][header]
						seq_out.write(record.format("fasta"))
						summary_out.write(header+"\t" + subset + "\n")
				outAlignment = msa.run_mafft(outdir +"/"+  subset + ".fa")
				outAlignmentFile = outdir +"/"+  subset + ".mafft.aln"
				outAlignmentWrite = SeqIO.write(outAlignment, outAlignmentFile, alignment_format)
			
		

		


		


		

		
		
		

	elif subcommand == "hmm_v_hmm":
		alignment_file1 = duomolog_args.aln1.name
		alignment_file2 = duomolog_args.aln2.name
		



