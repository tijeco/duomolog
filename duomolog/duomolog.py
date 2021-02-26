import sys
import argparse
from Bio import AlignIO
from Bio import SeqIO


from duomolog import blast



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
		
		parser.add_argument("-aln","-a", type=argparse.FileType('r'),required=True)
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
		

def notRedundant(sequences,aligned=True):
	seqs = {}
	for record in sequences:
		if record.seq not in seqs:
			if aligned:
				seq = str(record.seq).replace("-","")
			seqs[seq] = record.id
		else:
			seqs = {}
			break
	return seqs

def pairRedundant(known,unknown, aligned=True):
	nr_known = notRedundant(known)
	nr_unknown = notRedundant(unknown)
	
	if nr_known == False:
		sys.exit("input alignment contains redundant sequences, please fix and try again")
	if nr_unknown == False:
		sys.exit("input query contains redundant sequences, please fix and try again")
	
	known_seqs = set(nr_known.keys())
	unknown_seqs = set(nr_unknown.keys())
	shared_seqs = known_seqs & unknown_seqs
	shared_seqs_unknown_ids = [nr_unknown[i] for i in shared_seqs]
	return shared_seqs_unknown_ids
	
	
	
	

	# print("checking if ", known, "has ")
	# if aligned:
		







def main():
	cli = ParseCommands()

	duomolog_args = cli.args
	subcommand = cli.base_args.command 
	alignment_format = cli.base_args.format

	query_file = duomolog_args.query.name
	query = SeqIO.parse(query_file, "fasta")

	if subcommand == "blast_v_hmmer":
		alignment_file = duomolog_args.aln.name
		alignment = AlignIO.read(open(alignment_file), alignment_format)
		sharedIDs = pairRedundant(alignment,query)
		if sharedIDs:
			print("Error: the following sequences were found in the query that are already in the input alignment, please remove and try again")
			for h in sharedIDs:
				print(h)
			sys.exit()
		
		
		

	elif subcommand == "hmm_v_hmm":
		alignment_file1 = duomolog_args.aln1.name
		alignment_file2 = duomolog_args.aln2.name
		



