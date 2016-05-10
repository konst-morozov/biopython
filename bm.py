
from bm_preproc import BoyerMoore

import bisect
   
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

class Index(object):
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer
    
    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def naive_2mm(p, t):
	# up to 2 mismatches allowed
	occurrences = []
	for i in range(len(t) - len(p) + 1):  # loop over alignments
		match = 2
		for j in range(len(p)):  # loop over characters
			if t[i+j] != p[j]:  # compare characters
				match -= 1
				if match == -1:
					break

		if match >= 0:
			occurrences.append(i)  # all chars matched; record
	return occurrences


def boyer_moore_with_counts(p, p_bm, t):
	""" Do Boyer-Moore matching. p=pattern, t=text,
		p_bm=BoyerMoore object for p """
	i = 0
	occurrences = []
	num_character_comparisons = 0
	num_alignments = 0

	while i < len(t) - len(p) + 1:
		shift = 1
		mismatched = False
		for j in range(len(p)-1, -1, -1):
			num_character_comparisons += 1
			if p[j] != t[i+j]:
				skip_bc = p_bm.bad_character_rule(j, t[i+j])
				skip_gs = p_bm.good_suffix_rule(j)
				shift = max(shift, skip_bc, skip_gs)
				mismatched = True
				break
		num_alignments += 1
		if not mismatched:
			occurrences.append(i)
			skip_gs = p_bm.match_skip()
			shift = max(shift, skip_gs)
		i += shift
	return occurrences, num_alignments, num_character_comparisons

	
def naive(p, t):
	occurrences = []
	for i in range(len(t) - len(p) + 1):  # loop over alignments
		match = True
		for j in range(len(p)):  # loop over characters
			if t[i+j] != p[j]:  # compare characters
				match = False
				break
		if match:
			occurrences.append(i)  # all chars matched; record
	return occurrences
	
def naive_with_counts(p, t):
	occurrences = []
	num_character_comparisons = 0
	num_alignments = 0
	for i in range(len(t) - len(p) + 1):  # loop over alignments
		match = True
		for j in range(len(p)):  # loop over characters
			num_character_comparisons += 1
			if t[i+j] != p[j]:  # compare characters
				match = False
				break
		num_alignments += 1
		if match:
			occurrences.append(i)  # all chars matched; record
	return occurrences, num_alignments, num_character_comparisons

def readGenome(filename):
	genome = ''
	with open(filename, 'r') as f:
		for line in f:
			# ignore header line with genome information
			if not line[0] == '>':
				genome += line.rstrip()
	return genome	
	
p = 'word'
t = 'there would have been a time for such a word'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)

p = 'needle'
t = 'needle need noodle needle'
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)


p = 'word'
t = 'there would have been a time for such a word'
lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)


p = 'needle'
t = 'needle need noodle needle'
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)


genome = readGenome('chr1.GRCh38.excerpt.fasta')
#_,al,_ = naive_with_counts('GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG', genome)
#print ("Number of alignments 1", al)

#_,_,cc = naive_with_counts('GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG', genome)
#print ("Number of occurrences 2", cc)

#p_bm = BoyerMoore('GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG', 'ACGT')
#_,al,_ = boyer_moore_with_counts('GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG', p_bm, genome)
#print ("Number of alignments 3", al)

#Number of alignments 1 799954
#Number of occurrences 2 984143
#Number of alignments 3 127974

#oc = naive_2mm('GGCGCGGTGGCTCACGCCTGTAAT',genome);
#print ("Number of occurrences (2mm) 4", len(oc))
#19

'''
Index-assisted approximate matching. In practicals, we built a Python class called Index
implementing an ordered-list version of the k-mer index. The Index class is copied below.

We also implemented the pigeonhole principle using Boyer-Moore as our exact matching algorithm.

Implement the pigeonhole principle using Index to find exact matches for the partitions. Assume P always has length 24, and that we are looking for approximate matches with up to 2 mismatches (substitutions). We will use an 8-mer index.

Download the Python module for building a k-mer index.

https://d28rh4a8wq0iu5.cloudfront.net/ads1/code/kmer_index.py

Write a function that, given a length-24 pattern P and given an Index object built on 8-mers, finds all approximate occurrences of P within T with up to 2 mismatches. Insertions and deletions are not allowed. Don't consider any reverse complements.

How many times does the string GGCGCGGTGGCTCACGCCTGTAAT, which is derived from a human Alu sequence, occur with up to 2 substitutions in the excerpt of human chromosome 1? (Don't consider reverse complements here.)

Hint 1: Multiple index hits might direct you to the same match multiple times, but be careful not to count a match more than once.
'''


'''
Using the instructions given in Question 4, how many total index hits are there when searching for occurrences of GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt of human chromosome 1?

(Don't consider reverse complements.)

Hint: You should be able to use the boyer_moore function (or the slower naive function) to double-check your answer.
'''



KMER = 8

def mm_num(t, p):
	mm = 0
	for i in range(len(p)):
		if t[i] != p[i]:
			mm +=1
	return mm

	
# TODO the below should be refactored#	
	
hits = []
ind = Index(genome, KMER)
hit1 = ind.query('GGCGCGGT')
if len(hit1) > 0:
	for v in hit1:
		if mm_num(genome[v+KMER:v+KMER*2], 'GGCTCACG') + mm_num(genome[v+KMER*2:v+KMER*3], 'CCTGTAAT') < 3:
			hits.append(v);
			
print(hit1)

hit2 = ind.query('GGCTCACG')
if len(hit2) > 0:
	for v in hit2:
		if mm_num(genome[v+KMER:v+KMER*2], 'CCTGTAAT') + mm_num(genome[v-KMER:v], 'GGCGCGGT') < 3:
			hits.append(v-KMER);
print(hit2)
			
hit3 = ind.query('CCTGTAAT')
if len(hit3) > 0:
	for v in hit3:
		if mm_num(genome[v-KMER*2:v-KMER], 'GGCGCGGT') + mm_num(genome[v-KMER:v], 'GGCTCACG') < 3:
			hits.append(v-KMER*2);
print(hit3)		
			
print('Number of times string appear (8-mer index)', len(set(hits)))

#print("naive 1:", len(naive('GGCGCGGT',genome))," 2:", len(naive('GGCTCACG',genome)), " 3:", len(naive('CCTGTAAT',genome)) )

print("Overall index hits:", len(hit1) + len(hit2) + len(hit3) )


'''
Let's examine whether there is a benefit to using an index built using subsequences of T rather than substrings, as we discussed in the "Variations on k-mer indexes" video. We'll consider subsequences involving every N characters. For example, if we split ATATAT into two substring partitions, we would get partitions ATA (the first half) and TAT (second half). But if we split ATATAT into two subsequences by taking every other character, we would get AAA (first, third and fifth characters) and TTT (second, fourth and sixth).

Another way to visualize this is using numbers to show how each character of P is allocated to a partition. Splitting a length-6 pattern into two substrings could be represented as 111222, and splitting into two subsequences of every other character could be represented as 121212
The following class SubseqIndex is a more general implementation of Index that additionally handles subsequences. It only considers subsequences that take every Nth character:


For example, if we do:

ind = SubseqIndex('ATATAT', 3, 2)
print(ind.index)
we see:

[('AAA', 0), ('TTT', 1)]
And if we query this index:

p = 'TTATAT'
print(ind.query(p[0:]))
we see:

[]
because the subsequenceTAAis not in the index. But if we query with the second subsequence:

print(ind.query(p[1:]))
we see:

[1]
because the second subsequence TTT is in the index.

Write a function that, given a length-24 pattern P and given a SubseqIndex object built with k = 8 and ival = 3, finds all approximate occurrences of P within T with up to 2 mismatches.

When using this function, how many total index hits are there when searching for GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt of human chromosome 1? (Again, don't consider reverse complements.)

Hint: See this notebook for a few examples you can use to test your function.
'''


def getHits(t, p, hits, offset):
	occur = []
	#print (hits)
	if len(hits) > 0:
		for v in hits:  
			if (mm_num(t[v-offset: v+24-offset], p) < 3):
				occur.append(v-offset);
	return occur

def query_subseq(p, t, subseq_ind):
	hits = []
	hit1 = subseq_ind.query(p[0:])
	hit2 = subseq_ind.query(p[1:])
	hit3 = subseq_ind.query(p[2:])
	
	hits.extend(getHits(t, p, hit1, 0))
	hits.extend(getHits(t, p, hit2, 1))
	hits.extend(getHits(t, p, hit3, 2))
	
	return set(hits), len(hit1)+len(hit2)+len(hit3)
	
	
t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
subseq_ind = SubseqIndex(t, 8, 3)

occurrences, num_index_hits = query_subseq(p, t, subseq_ind)
print("tomorrow hits:", occurrences)
#[0, 14]

print("tomorrow index:",num_index_hits)
#6



# King John by William Shakespeare
# !wget http://www.gutenberg.org/ebooks/1110.txt.utf-8

#bitsadmin  /transfer mydownloadjob  /download  /priority normal  http://www.gutenberg.org/ebooks/1110.txt.utf-8  C:\Users\morozko\Downloads\1110.txt.utf-8

#powershell -command "& { iwr http://www.gutenberg.org/cache/epub/1110/pg1110.txt  -OutFile 1110.txt.utf-8 }"

#import wget
#filename = wget.download('http://www.gutenberg.org/ebooks/1110.txt.utf-8')


t = open('pg1110.txt').read()
p = 'English measure backward'

subseq_ind = SubseqIndex(t, 8, 3)

occurrences, num_index_hits = query_subseq(p, t, subseq_ind)

print(occurrences)
#[135249] OR 132185 :) file format UTF8 or whatever
print(num_index_hits)
#3

subseq_ind = SubseqIndex(genome, 8, 3)
_, num_index_hits = query_subseq("GGCGCGGTGGCTCACGCCTGTAAT", genome, subseq_ind)
print('number of hits: ', num_index_hits)

