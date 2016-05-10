def overlap(a, b, min_length=3):
	""" Return length of longest suffix of 'a' matching
		a prefix of 'b' that is at least 'min_length'
		characters long.  If no such overlap exists,
		return 0. """
	start = 0  # start all the way at the left
	while True:
		start = a.find(b[:min_length], start)  # look for b's suffx in a
		if start == -1:  # no more occurrences to right
			return 0
		# found occurrence; check for full suffix/prefix match
		if b.startswith(a[start:]):
			return len(a)-start
		start += 1  # move just past previous match

import itertools

def scs(ss):
	""" Returns shortest common superstring of given
		strings, which must be the same length """
	shortest_sup = ['AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA']
	for ssperm in itertools.permutations(ss):
		sup = ssperm[0]  # superstring starts as first string
		for i in range(len(ss)-1):
			# overlap adjacent strings A and B in the permutation
			olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
			# add non-overlapping portion of B to superstring
			sup += ssperm[i+1][olen:]
		if len(sup) < len(shortest_sup[0]):
			shortest_sup = [sup]  # found shorter superstring
		elif len(sup) == len(shortest_sup[0]):
			shortest_sup.append(sup);
	return shortest_sup  # return shortest

'''
It's possible for there to be multiple different shortest common superstrings for the same set of input strings. Consider the input strings ABC, BCA, CAB. One shortest common superstring is ABCAB but another is BCABC and another is CABCA.

What is the length of the shortest common superstring of the following strings?
CCT, CTT, TGC, TGG, GAT, ATT
'''
	

ss= ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'];

import itertools
list = []
'''
for sss in itertools.permutations(ss):
#	print(sss)
	list.extend(scs(sss));
'''	
print ('SS', set(list))
#11

print ('Len ', len(set(list)))
#4

'''
How many different shortest common superstrings are there for the input strings given in the previous question?

Hint 1: You can modify the scs function to keep track of this.

Hint 2: You can look at these examples to double-check that your modified scs is working as expected.

'''

strings = ['ABC', 'BCA', 'CAB']
# Returns list of all superstrings that are tied for shorest
print(scs(strings))
#['ABCAB', 'BCABC', 'CABCA']


strings = ['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']
print(scs(strings))
'''
['AATAGATCGTGC',
 'AATAGATGCTCG',
 'AATAGTCGATGC',
 'AATCGATAGTGC',
 'AATGCTCGATAG',
 'TCGAATAGATGC',
 'TCGATAGAATGC',
 'TCGATGCAATAG',
 'TGCAATAGATCG',
 'TGCAATCGATAG']
'''


'''
Download this FASTQ file containing synthetic sequencing reads from a mystery virus:

https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ads1_week4_reads.fq

All the reads are the same length (100 bases) and are exact copies of substrings from the forward strand of the virus genome. You don't have to worry about sequencing errors, ploidy, or reads coming from the reverse strand.

Assemble these reads using one of the approaches discussed, such as greedy shortest common superstring. Since there are many reads, you might consider ways to make the algorithm faster, such as the one discussed in the programming assignment in the previous module.

How many As are there in the full, assembled genome?

Hint: the virus genome you are assembling is exactly 15,894 bases long
'''


'''
How many Ts are there in the full, assembled genome from the previous question?
'''

#ROOT OF N - some precision 

import math
A = 1244444.0
B = A
for i in range(30):
	A = math.sqrt(A)
	
A -=1
A = A/9
A +=1

for i in range(30):
	A = A**2

print(A*A*A*A*A*A*A*A*A - B)

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def pick_maximal_overlap(reads, k, pairs):
	reada, readb = None, None
	best_olen = 0
	for a,b in itertools.permutations(reads, 2):
		olen = overlap (a, b, min_length = k)
		if olen > best_olen:
			reada, readb = a, b
			best_olen = olen
	return reada, readb, best_olen
	
def greedy_scs(reads, k):
	read_a, read_b, olen = pick_maximal_overlap(reads, k)
	while olen > 0:
		reads.remove(read_a)
		reads.remove(read_b)
		reads.append(read_a+read_b[olen:])
		print(len(reads))
		read_a, read_b, olen = pick_maximal_overlap(reads, k)
	return ''.join(reads)

	
def overlap_all_pairs(reads, kmer):
	pairs = []
	dict = {}
	for read in reads:
		for i in range(len(read) - kmer + 1):
			k = read[i: i + kmer]
			if k in dict: 
				dict[k].append(read)
			else:
				dict[k] = [read]
	#print ('D', dict)
	
	for read in reads:
		suf = read[-kmer:]
#		print(suf)
	
		for elem in dict[suf]:
			if elem != read:
				if overlap(read, elem, kmer) != 0:
					pairs.append((read, elem))
	
	return set(pairs);
	

def pick_maximal_overlap_with_pairs(reads, k):
	reada, readb = None, None
	best_olen = 0
	pairs = overlap_all_pairs(reads, k)
	
	for pair in pairs:
		olen = overlap (pair[0], pair[1], min_length = k)
		if olen > best_olen:
			reada, readb = pair[0], pair[1]
			best_olen = olen
	return reada, readb, best_olen

	
def greedy_scs_with_pairs(reads, k):
	read_a, read_b, olen = pick_maximal_overlap_with_pairs(reads, k)
	while olen > 0:
		reads.remove(read_a)
		reads.remove(read_b)
		reads.append(read_a+read_b[olen:])
		print(len(reads))
		read_a, read_b, olen = pick_maximal_overlap_with_pairs(reads, k)
	return ''.join(reads)

	
#print(greedy_scs(['ABCD', 'CDBC', 'BCDA'], 1))


def overlap_all_pairs(reads, kmer):
	pairs = []
	dict = {}
	for read in reads:
		for i in range(len(read) - kmer + 1):
			k = read[i: i + kmer]
			if k in dict: 
				dict[k].append(read)
			else:
				dict[k] = [read]
	#print ('D', dict)
	
	for read in reads:
		suf = read[-kmer:]
#		print(suf)
	
		for elem in dict[suf]:
			if elem != read:
				if overlap(read, elem, kmer) != 0:
					pairs.append((read, elem))
	
	return set(pairs);
	
def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities
	
	
	
reads,_ = readFastq('ads1_week4_reads.fq')

str = greedy_scs_with_pairs(reads, 10)

print('LEN: ', len(str), '. A: ', str.count ('A'),  '. T: '  , str.count ('T'))

#LEN:  15894 . A:  4633 . T:  3723