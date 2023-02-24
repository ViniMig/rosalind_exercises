#########################################################################################################
# This python file contains all the function definitions used mainly for solving Rosalind exercises.    #
# A few functions present are actually answers to other courses and a comment will inform so.           #
# From time to time an effort is made towards refactoring and improving my coding skills.               #
# It is also being updated  with new exercises solutions.                                               #
#########################################################################################################

from itertools import permutations, product
from collections import defaultdict, Counter
from typing import List
from scipy.special import binom
import math
import re
import pandas as pd
import numpy as np
import requests
import bisect
import networkx as nx

def readFasta(f: str) -> dict:
    """
    Reads a fasta file and
    Returns a dictionary with ID as key and sequence as value.
    """
    fastaFile = open(f, 'r')
    fastaLines = fastaFile.readlines()
    seqID: str = ""
    seq: str = ""
    fastaSeqs: dict = {}
    for fasta in fastaLines:
        if fasta.startswith(">"):
            fastaSeqs[seqID] = seq
            seqID = fasta.rstrip().removeprefix(">")
            seq = ""
        else:
            seq += fasta.rstrip()
    fastaSeqs[seqID] = seq
    fastaSeqs.pop("")
    fastaFile.close()
    
    return fastaSeqs

def count_nuclts(s:str) -> dict:
    """
    Returns the number of nucleotides in a given sequence.
    """
    return Counter(s)

def transcribe_dna(t: str) -> str:
    """
    Returns the sequence 't' transcribed into RNA.
    """
    return(t.replace('T', 'U'))

def rev_complement(s: str) -> str:
    """
    Returns the Reverse Complement of a Given sequence 's'.
    """
    s  = s[::-1]
    s = s.replace("A","%temp%").replace("T", "A").replace("%temp%", "T").replace("C", "%temp%").replace("G","C").replace("%temp%","G")
    return(s)

def rabbit_seq(n: int, k:int) -> int:
    """
    Returns the Rabbit population 'n' after 'k' generations.
    """
    if n <= 2:
        return 1
    else:
        return ((k * rabbit_seq(n-2,k)) + rabbit_seq(n-1,k))

def highest_gc_content(fp: str) -> tuple:
    """
    Returns the Sequence ID and the GC content in percentage of the
    Given 'fp' fasta file.
    """
    highGC: float = 0
    highID: str = ""

    fasta_dic: dict = readFasta(fp)
    for k, val in fasta_dic.items():
        count_curr_nuctds: dict = count_nuclts(val)
        curr_gc_cont: float = ( (count_curr_nuctds["G"] + count_curr_nuctds["C"]) / len(val) ) * 100
        if curr_gc_cont > highGC:
            highGC = curr_gc_cont
            highID = k
    return highID, highGC

def hamming_dist(s: str, t: str) -> int:
    """
    Returns the Hamming distance between two
    Given sequences 's' and 't'.
    """
    hamdist: int = 0

    for (base_a, base_b) in zip(s,t):
        if base_a != base_b:
            hamdist += 1
    
    return hamdist

def dominant_allele(k: int, m: int, n: int) -> float:
    """
    Returns the probability that 2 randomly selected organisms will produce an individual possessing a dominant allele.
    Population is made of 'k' homozygous dominant + 'm' heterozygous + 'n' homozygous recessive.
    """
    pop: int = k + m + n
    prob_no_dominant: float = ( ((n**2) - n) + (((m**2) - m) * 0.25) + (m * n) ) / ((pop**2) - pop)

    return (1 - prob_no_dominant)

################################################
# Create dictionary with rna codons
# and the respective protein translated by each.
# File of type:
# codon   aminoacid
# AAA   A
# CCC   B
# TTT   C
# GGG   D
prot_dic: dict = {}
with open('codon_tab.txt', 'r') as prot_file:
    codons = prot_file.readlines()
    for codon in codons:
        key_val: list = codon.rstrip().split(" ")
        prot_dic[key_val[0]] = key_val[1]
################################################

def translate_rna(s: str) -> str:
    """
    For a given rna sequence 's'
    Returns the translated sequence until the stop codon is found.
    It does not take into account the Start Codon
    """
    prot_seq: str = ""
    mrna_reads = re.findall("[ACTG]{3}",s)
    for mrna_codon in mrna_reads:
        if prot_dic[mrna_codon] == "Stop":
            break
        else:
            prot_seq += prot_dic[mrna_codon]

    return prot_seq

def find_motif(s: str, t: str) -> list:
    """
    Finds the motif 't' in the sequence 's'.
    Prints the locations.
    """
    locales = re.finditer(r'(?=('+t+'))', s)
    return [local.span()[0] + 1 for local in locales]

def consensus_profile(fp: str) -> tuple:
    """
    Returns a consensus sequence and a profile matrix for a
    Given set of sequences in fasta format.
    """
    fastaDic: dict = readFasta(fp)
    num_seqs = len(fastaDic.values())
    seqs_matrix = np.array(list(list(fastaDic.values())[0]))
    for value in fastaDic.values():
        seqs_matrix = np.vstack((seqs_matrix, np.array(list(value))))
    seqs_matrix = seqs_matrix.astype(np.matrix)
    seqs_matrix = np.delete(seqs_matrix, obj = 0, axis=0)
    seqs_df = pd.DataFrame(seqs_matrix) #turn into pandas dataframe for easier manipulation
    
    #determine consensus
    consensus: str = ''.join(seqs_df[column].mode()[0] for column in seqs_df.columns)

    aa_counts = seqs_df.apply(pd.value_counts, axis = 0).fillna(0)
    aa_counts = np.matrix(aa_counts) #ACGT -> row for each nuc
    aa_entropy = aa_counts / num_seqs
    aa_entropy = np.multiply(aa_entropy, np.where(aa_entropy != 0, np.log2(aa_entropy), 0))
    aa_entropy = np.sum(aa_entropy.sum(axis=0) * -1)
    np.savetxt("counts_matrix.txt", aa_counts, fmt='%d',)
    return(consensus, aa_counts, aa_entropy)

# Initial dictionary for mortal_rabbit function.
# Month 1 and Month 2 assume the rabbits are born and need
# to mature for 1 month. The value for 0 is for the fact
# that the function is recursive.
rabbit_pops = {
        0: 1,
        1: 1,
        2: 1
        }

def mortal_rabbit(n: int, m:int) -> int:
    """
    Returns the total pairs of rabbits that will remain after
    'n' months if all rabbits live for 'm' months. Using dynamic programing, so we don't need to
    recalculate every time.
    """
    if n not in rabbit_pops.keys():
        if n <= m:
            rabbit_pops[n] = mortal_rabbit(n-1, m) + mortal_rabbit(n-2, m)
        else:
            rabbit_pops[n] = mortal_rabbit(n-1, m) + mortal_rabbit(n-2, m) - mortal_rabbit(n - (m + 1), m)
    
    return rabbit_pops[n]

def over_graphs(fasta: str) -> list:
    """
    Returns the adjacency list for the
    Given 'fasta' sequences.
    """
    fastas: dict = readFasta(fasta)
    adjList: list = []

    for currentID, currentSeq in fastas.items():
        currentList: list = [currentID + " " + k for k,v in fastas.items() if k != currentID and v.startswith(currentSeq[-3::])]
        adjList.extend(currentList)

    return(adjList)

def expected_dominant(nums: str) -> float:
    """
    Returns the expected number of offspring with dominant phenotypes in the next generation
    Given six integers corresponding to the number of couples posessing each genotype pairing: \n
    AA-AA    AA-Aa    AA-aa    Aa-Aa    Aa-aa    aa-aa \n
    Every couple is assumed to produce two offspring
    """
    nums_list: list = [int(num) for num in nums.split(" ")]
    expected: float = ( (2 * nums_list[0] * 1) + (2 * nums_list[1] * 1) + (2 * nums_list[2] * 1) + (2 * nums_list[3] * 0.75) + (2 * nums_list[4] * 0.5) )
    return expected

def longest_common(fastaPath: str) -> str:
    """
    Retuns the longest common motif shared by all the sequences in the
    Given 'fastaPath' file.
    """
    fasta_dict: dict = readFasta(fastaPath)

    shortest_seq: str = sorted(list(fasta_dict.values()), key=len)[0] #the longest motif will be at most the size of the smallest sequence.

    # get all possibilities of motifs, slicing the string in all the possible ways
    # sort from longest to shortest and stop as soon as we find a common motif to all the other sequences, this will be the longest possible.
    possible_motifs: list = [shortest_seq[i: j] for i in range(len(shortest_seq)) for j in range(i + 1, len(shortest_seq) + 1)]
    possible_motifs = sorted(set(possible_motifs), key= len, reverse=True) #remove duplicates
    for motif in possible_motifs:
        if all(motif in other_seq for other_seq in fasta_dict.values()):
            return motif
    
    return "No shared motifs."

#########################################################################################
## clump exercises are for the Coursera Bioinformatics course.
# def get_clump_freqs(clump: str, seq_size: int):
#     freqs_dict = dict()
#     possible_seqs = [clump[i:i+seq_size] for i in range(0, len(clump) - seq_size)]
#     for seq in possible_seqs:
#         freqs_dict[seq] = freqs_dict.get(seq, 0) + 1
#     return freqs_dict

# def find_clump(g_nome: str, k: int, L: int, t: int):
#     clumpWindows = [g_nome[i:i+L] for i in range(0,len(g_nome) - L)]
#     clump_L_t = []
#     for clump in clumpWindows:
#         counts = get_clump_freqs(clump, k)
#         clump_L_t += [w for (w,v) in counts.items() if v >= t and w not in clump_L_t]

#     return clump_L_t
## end of clump functions
#########################################################################################

def independent_alleles(k: int, N:int) -> float:
    """
    Returns Prob of at least 'N' offspring are AaBb from 'k'th generation.
    """
    pop: int = 2**k
    prob: float = 0

    for n in range(N, pop + 1):
        prob += (math.factorial(pop) / (math.factorial(n) * math.factorial(pop - n))) * (0.25** n) * (0.75**(
                pop - n))
    return prob

def independent_alleles_alt(k: int, N:int) -> float:
    """
    Returns Prob of at least 'N' offspring are AaBb from 'k'th generation.
    """
    pop: int = 2**k
    prob: float = 0

    for n in range(N, pop + 1):
        prob += (binom(pop, n)) * (0.25** n) * (0.75**(
                pop - n))
    return prob

#motifs

def retrieve_uniprot_fasta(uniprot_id: str) -> str:
    """
    Retrieves the uniprot fasta using the requests module.
    Returns a string with the fasta sequence of the protein.
    """
    uniprot_id = uniprot_id.split('_', 1)[0]

    url: str = f'https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta'
    get_fasta: str = requests.get(url).text
    get_fasta = get_fasta.split('\n', 1)[1]
    get_fasta = ''.join(get_fasta.rstrip().split('\n'))
    
    return get_fasta

def retrieve_fastas(ids_file: str) -> dict:
    """
    Reads the input file with the uniprot IDs and for each
    starts a request on uniprot.
    Returns a dictionary where uniprot ids are the keys and
    fasta text are the values.
    """
    with open(ids_file, "r") as f:
        uniprot_data: dict = {line.rstrip(): retrieve_uniprot_fasta(line.rstrip()) for line in f.readlines()}
    
    return uniprot_data

def motif_location(fasta_dict: dict) -> None:
    """
    Receives a dictionary of ids and their fastas and
    Prints the ids and their N-glycosylation 
    motifs if they exist.
    """
    n_glyco_expr: str = "N[^P][ST][^P]"
    for k,v in fasta_dict.items():
        print(k)
        find_motif(v, n_glyco_expr)
        
#mRNA from prota
##############################################################################
# Create reverse dict of the codons list
# File of type:
# codon   aminoacid
# AAA   A
# CCC   B
# TTT   C
# GGG   D
protein_seqs = dict()
with open("codon_tab.txt", "r") as f:
    for codons in f.readlines():
        codon_prot: list = codons.rstrip().split()
        protein_seqs[codon_prot[1]] = protein_seqs.get(codon_prot[1], list())
        protein_seqs[codon_prot[1]].append(codon_prot[0])
##############################################################################

def count_possible_rnas(protein_seq: str) -> int:
    """
    For the given input Protein Sequence
    Returns the total number of possible RNA strings modulo 1,000,000
    """
    total_seqs: int = len(protein_seqs["Stop"]) % 1000000 #account for the possible stop codons

    for aa in protein_seq:
        total_seqs = (total_seqs * len(protein_seqs[aa]) % 1000000) %1000000
    
    return total_seqs

#orfs

def translate_rna_with_stop(s: str) -> str:
    """
    For a Given rna sequence
    Returns the translated sequence.
    It does not take into account the Start Codon, and Stop is included.
    """
    mrna_reads = re.findall("[ACGU]{3}",s)
    prot_seq: str = ''.join(prot_dic[rna_codon] for rna_codon in mrna_reads)

    return prot_seq

def possible_seqs(rna_seq: str) -> list:
    """
    Returns a protein sequence for the given RNA seq.
    Account for start codon M.
    """
    prot_seq_list = []
    
    prot_seq: list = translate_rna_with_stop(rna_seq).split("Stop")
    prot_seq.pop() #drop last sequence, it will either be empty (meaning it ends with a Stop) or some sequence without Stop codon
    for protein_seq in prot_seq:
        if "M" in protein_seq:
            pos_M = re.finditer('M', protein_seq)
            for position_M in pos_M:
                idx = position_M.span()[0]
                prot_seq_list.append(protein_seq[idx:])

    return prot_seq_list

def or_frames(fasta_file: str) -> set:
    """
    Returns all the possible protein sequences for the given fasta sequence.
    """
    fastaF: dict = readFasta(fasta_file)

    current_seq: str = list(fastaF.values())[0]
    reverse_compl: str = rev_complement(current_seq)
    transc_seq: str = transcribe_dna(current_seq)
    transc_rev_seq: str = transcribe_dna(reverse_compl)

    possible_sequences = []
    for i in range(3):
        possible_sequences.extend(possible_seqs(transc_seq[i:]))
        possible_sequences.extend(possible_seqs(transc_rev_seq[i:]))

    return set(possible_sequences)

#permutations
def permutations_rosalind(n: int) -> tuple:
    list_perms = list(permutations(range(1, n+1)))
    with open("perms.txt", "w") as f:
        f.write(str(len(list_perms)))
        f.write("\n")
        for seq in list_perms:
            f.write(' '.join(str(e) for e in seq))
            f.write("\n")

#permutations_rosalind(5)

#Protein Mass

###########################################################################################################
# Creqate dictionary of aminoacid masses from existing file.
# File of type:
# aminoacid   mass
# A   123.45678
# B   12.34567
with open("monoisotopic_mass.txt", "r") as f:
    mass_dict: dict = {line.rstrip().split()[0]: float(line.rstrip().split()[1]) for line in f.readlines()}
###########################################################################################################

def protein_mass(protein_sequence: str) -> float:
    """
    Returns the sum of the individual monoisotopic masses of each aminoacid plus one water for the input sequence.
    """
    return sum(mass_dict[aa] for aa in protein_sequence)

#Restriction sites

def test_palindrome(seq_a: str, seq_b: str, seq_size: int) -> list[tuple[int]]:
    """
    Returns list of pairs position and seq_size of the existing palindromes of given 2 sequences and seq_size.
    """
    palindromes = []
    for i in range(len(seq_a) - seq_size +1):
        if i == 0:
            if seq_a[:seq_size] == seq_b[(seq_size) - 1::-1]:
                palindromes.append((i + 1, seq_size))
        elif seq_a[i:i+seq_size] == seq_b[(i + seq_size) - 1: (i - 1):-1]:
            palindromes.append((i + 1, seq_size))
    return palindromes


def reverse_palindromes(fasta: str) -> None:
    """
    Gets the path to the fasta file containing the target sequence and
    Writes to the file 'palindrome.txt' the locations and length of the possible palindromes.
    """
    dna_seq = list(readFasta(fasta).values())[0] #only one fasta sequence given
    dna_compl: str = rev_complement(dna_seq)[-1::-1]
    pal_list = []
    for i in range(4, 13):
        pal_list.extend(test_palindrome(dna_seq, dna_compl, i))
    with open("palindromes.txt", "w") as res:
        for a,b in pal_list:
            res.write(f"{a}\t{b}\n")

#RNA Splicing

def clean_introns(dna_seq: str, introns: list[str]) -> str:
    """
    For given sequence string and intron list removes the introns and
    Returns the coding region of the sequence.
    """
    for seq in introns:
        dna_seq = dna_seq.replace(seq, '')
    return dna_seq

#Enumerate k-mers
# TODO:maybe write the output to a file instead.
def combine_and_sort(symbl: str, n: int) -> None:
    """
    Reads a sequence of symbols, splits into a list and combines them into sequences of length n.
    Prints the sequences sorted alphabetically
    """
    symbol_list = symbl.split()
    possible_combs = sorted(list(product(symbol_list, repeat = n)))
    for comb in possible_combs:
        print(''.join(comb))

#Longest ordered subsequence
# Implementing solution found on: https://leetcode.com/problems/longest-increasing-subsequence/solutions/1326308/c-python-dp-binary-search-bit-segment-tree-solutions-picture-explain-o-nlogn/

def pathOfLIS(nums: List[int]):
    sub = []
    subIndex = []  # Store index instead of value for tracing path purpose
    trace = [-1] * len(nums)  # trace[i] point to the index of previous number in LIS
    for i, x in enumerate(nums):
        if len(sub) == 0 or sub[-1] < x:
            if subIndex:
                trace[i] = subIndex[-1]
            sub.append(x)
            subIndex.append(i)
        else:
            idx = bisect.bisect_left(sub, x)  # Find the index of the smallest number >= x, replace that number with x
            if idx > 0:
                trace[i] = subIndex[idx - 1]
            sub[idx] = x
            subIndex[idx] = i

    path = []
    t = subIndex[-1]
    while t >= 0:
        path.append(nums[t])
        t = trace[t]
    return path[::-1]

def findLDS(nums: List[int]):

    # base case
    if not nums:
        return

    # `LDS[i]` stores the longest decreasing subsequence of sublist
    # `nums[0…i]` that ends with `nums[i]`
    LDS = [[] for _ in range(len(nums))]

    # `LDS[0]` denotes longest decreasing subsequence ending at `nums[0]`
    LDS[0].append(nums[0])

    # start from the second element in the list
    for i in range(1, len(nums)):
        # do for each element in sublist `nums[0…i-1]`
        for j in range(i):
            # find longest decreasing subsequence that ends with `nums[j]`
            # where `nums[j]` is more than the current element `nums[i]`
            if nums[j] > nums[i] and len(LDS[j]) > len(LDS[i]):
                LDS[i] = LDS[j].copy()

        # include `nums[i]` in `LDS[i]`
        LDS[i].append(nums[i])

    # `j` will contain an index of LDS
    j = 0
    for i in range(len(nums)):
        if len(LDS[j]) < len(LDS[i]):
            j = i

    # print LDS
    return LDS[j]

#Genome Assembly

SeqList = List[str]

# Implementing solution from: https://www.youtube.com/watch?v=uS6ca7yeVb0

def overlap(a: str, b: str, min_length: int = 3) -> int:
    """
    Returns length of longest suffix of 'a' matching
    a prefix of 'b' that is at least 'min_length'
    characters long. If no such overlap exists,
    returns 0.
    """
    start: int = 0 #start all the way at the left
    while True:
        start = a.find(b[:min_length], start) # look for b's suffix in 'a' (returns the lowest index in a)
        if start == -1: #no occurences
            return 0
        #found occurence, check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1 #move past previous match

def max_overlap(list_s: SeqList, k: int = 3) -> tuple:
    """
    Determines the maximum overlapping 2 strings from the list with overlap size k and
    Returns the list updated with the maximum overlap
    """
    max_size: int = 0 #max size of overlap
    seq_1, seq_2 = None, None

    for a,b in permutations(list_s,2):
        max_over: int = overlap(a, b, min_length = k)
        if max_over > max_size:
            seq_1, seq_2 = a, b
            max_size = max_over
    
    return seq_1, seq_2, max_size

def gen_assembly(list_seqs: SeqList, k:int =3) -> str:
    """
    Receives a List of sequences and
    Returns the smallest superstring
    """
    read_a, read_b, olen = max_overlap(list_seqs, k)

    while olen > 0:
        list_seqs.remove(read_a)
        list_seqs.remove(read_b)
        list_seqs.append(read_a + read_b[olen:])
        read_a, read_b, olen = max_overlap(list_seqs, k)
    
    return ''.join(list_seqs)