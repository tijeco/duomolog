# duomolog
A method to identify the best set of homologous sequences from two homology searching approaches.

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
