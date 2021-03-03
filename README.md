# duomolog
A method to identify the best set of homologous sequences from two homology searching approaches.


# Usage
```
usage: duomolog <command> [<args>] <queryPepFasta>
                        The most commonly used duomolog commands are:
                                blast_v_hmmer
                                hmm_v_hmm

usage: duomolog blast_v_hmmer [<args>] 
Runs blast and hmmer

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Verified homologous input protein sequence fasta file
  --query QUERY, -q QUERY
                        FASTA formatted file containing database of peptides to be searched
  --intersect_only [INTERSECT_ONLY]
                        Only write hits from both approaches
  --blastout BLASTOUT   Tab delimited output file from BLAST
  --hmm HMM             HMM output file from hmmbuild
  --outFile OUTFILE, -o OUTFILE
                        Name of file to write output to
```
# Sub-commands
`blast_v_hmmer`
- compares hits of BLAST results vs HMMER

`hmm_v_hmm`
- Compares results from HMMER using two different Hidden Markov Models


# Duo comparison

The `Duo` object in `duo.py` takes in two `set()` objects of strings and returns the following comparisons (given a `left` set and a `right` set) with [python set operations](https://docs.python.org/3/library/stdtypes.html#set).
## Left
![union](figures/left.svg)
## Right
![union](figures/right.svg)

# Exclusive
## Left only
![union](figures/left_not_right.svg)
`left - right`
## Right only
![union](figures/right_not_left.svg)
`right - left`
# Intersection
![union](figures/left_intersect_right.svg)
`left & right`
# Union
![union](figures/left_union_right.svg)
`left | right`


# Pick the best subset

<!-- To be conservative, for now I am only focusing on the intersection. In the future, I would like to construct a multiple sequence alignment from the query  hits and the verified input, and compare the alignments and pick the subset that has the "best alignment" (less gappy, ) -->