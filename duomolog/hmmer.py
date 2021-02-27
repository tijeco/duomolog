import pyhmmer
import shlex, subprocess
import os

def run_hmmer(alnFile,queryFile):
	devnull = open(os.devnull, 'w')
	hmmFile = alnFile+".hmm"
	call_list = ''.join(['hmmbuild ',hmmFile,' ', alnFile])   
	commands = shlex.split(call_list)  
	subprocess.Popen(commands, stdin=subprocess.PIPE,           
        	stderr=subprocess.PIPE,stdout=devnull).communicate() 
	
	
	with pyhmmer.plan7.HMMFile(hmmFile) as h:
		hmm = next(h)
	with pyhmmer.easel.SequenceFile(queryFile) as seq_file:
		sequences = [ seq.digitize(hmm.alphabet) for seq in seq_file ]
	pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
	hits = pipeline.search_hmm(hmm, sequences) # Has lots of goodies!
	return hits
	
	