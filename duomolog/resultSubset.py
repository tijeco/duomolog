import sys



def seqSubSet(seq,subset):
	seqDict = {}
	for subsetName in subset:
		seqDict[subsetName] = {}
		for header in subset[subsetName]:
			seqDict[subsetName][header] = seq[header]
	return seqDict

		
			


def writeOut(outdir,duo_subset,querySeq,intersect_only):
	if intersect_only:
		duo_subset.intersectOnly()
	else:
		duo_subset.dropRedundant()
	duo_subset.dropEmpty()
	if bool(duo_subset):
		blast_hmmer_subsetSeqs = seqSubSet(querySeq,duo_subset.subsets)
		with open(outdir +"/duomolog_results.txt", "w") as summary_out:
			for subset in blast_hmmer_subsetSeqs:
				with open(outdir +"/"+  subset + ".fa","w") as seq_out:
					for header in blast_hmmer_subsetSeqs[subset]:
						record = blast_hmmer_subsetSeqs[subset][header]
						seq_out.write(record.format("fasta"))
						summary_out.write(header+"\t" + subset + "\n")
	else:
		sys.exit("Error: nothing to write")
	

	
		
	
		
	
	# if intersect_only:
		

	# 	# print(len(blast_hmmer_subsetSeqs["blast_intersect_hmmer"]))
	# 	if len(blast_hmmer_subsetSeqs["blast_intersect_hmmer"]) >0:
	# 		with open(outdir +"/duomolog_results.txt", "w") as summary_out:
	# 			with open(outdir +"/blast_intersect_hmmer.fa","w") as seq_out:
	# 				# for input_header in inputSeq:
	# 				# 	record = inputSeq[input_header]
	# 				# 	seq_out.write(record.format("fasta"))
	# 				for header in blast_hmmer_subsetSeqs["blast_intersect_hmmer"]:
	# 					record = blast_hmmer_subsetSeqs["blast_intersect_hmmer"][header]
	# 					seq_out.write(record.format("fasta"))
	# 					summary_out.write(header+"\tblast_intersect_hmmer\n")
	# 			outAlignment = msa.run_mafft(outdir +"/blast_intersect_hmmer.fa")
	# 			outAlignmentFile = outdir +"/blast_intersect_hmmer.mafft.aln"
	# 			outAlignmentWrite = SeqIO.write(outAlignment, outAlignmentFile, alignment_format)
	# else:
	# 	blast_hmmer_subset.dropEmpty()
	# 	blast_hmmer_subset.dropRedundant()
	# 	blast_hmmer_subsetSeqs = seqSubSet(querySeq,duo_subset)
	# 	with open(outdir +"/duomolog_results.txt", "w") as summary_out:
	# 		for subset in blast_hmmer_subsetSeqs:
	# 			with open(outdir +"/"+  subset + ".fa","w") as seq_out:
	# 				# for input_header in inputSeq:
	# 				# 	record = inputSeq[input_header]
	# 				# 	seq_out.write(record.format("fasta"))
	# 				for header in blast_hmmer_subsetSeqs[subset]:
	# 					record = blast_hmmer_subsetSeqs[subset][header]
	# 					seq_out.write(record.format("fasta"))
	# 					summary_out.write(header+"\t" + subset + "\n")
	# 			outAlignment = msa.run_mafft(outdir +"/"+  subset + ".fa")
	# 			outAlignmentFile = outdir +"/"+  subset + ".mafft.aln"
	# 			outAlignmentWrite = SeqIO.write(outAlignment, outAlignmentFile, alignment_format)






