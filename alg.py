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
		
def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t
	
	
def naive_with_rc(p, t):
	occurences = naive(p, t)
	complement = reverseComplement(p)
	if complement != p:
		occurences += naive(complement, t)
	return occurences
	

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def readGenomeLines(filename):
	genome = {}
	id = ''
	with open(filename, 'r') as f:
		for line in f:
            # ignore header line with genome information
			if line[0] == '>':
				id = line[1:].rstrip()
			elif len(line) > 0:
				if id in genome.keys():
					genome[id]+= line.rstrip()
				else: 
					genome[id] = line.rstrip()
	return genome
	

	
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

	

#print(reverseComplement('ACCAACAAC'))
#print(readGenome('lambda_virus.fa')[:7])

print(len(naive_with_rc("AGGT", readGenome('lambda_virus.fa'))))
print(len(naive_with_rc("TTAA", readGenome('lambda_virus.fa'))))

print(min(naive_with_rc("ACTAAGT", readGenome('lambda_virus.fa'))))
print(min(naive_with_rc("AGTCGA", readGenome('lambda_virus.fa'))))

print(len(naive_2mm("TTCAAGCC", readGenome('lambda_virus.fa'))))
print(min(naive_2mm("AGGAGGTT", readGenome('lambda_virus.fa'))))

p = 'CCC'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
occurrences = naive_with_rc(p, t)
print(occurrences)

p = 'CGCG'
t = ten_as + 'CGCG' + ten_as + 'CGCG' + ten_as
occurrences = naive_with_rc(p, t)
print(occurrences)

occurrences = naive_with_rc('ATTA', readGenome('phix.fa'))
print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))


p = 'CTGT'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as
occurrences = naive_2mm(p, t)
print(occurrences)

occurrences =  naive_2mm('GATTACA', readGenome('phix.fa'))
print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))


reads, quals = readFastq('ERR037900_1.first1000.fastq')
print(len(reads))
print(len(quals))


	

import sys

print(sys.version)

import random
def create_dna(n, alphabet='acgt'):
    return ''.join([random.choice(alphabet) for i in range(n)])
print(create_dna(2))


from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
my_seq = Seq("TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG")
#print(my_seq.translate())

#results_handle = NCBIWWW.qblast('blastn','nr',"TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG")
#blast_records = NCBIXML.parse(results_handle)
#blast_record = next(blast_records)

#for alignment in blast_record.alignments:
#    print('sequence:', alignment.title)


#genome = readGenomeLines('dna.example.fasta')
genome = readGenomeLines('dna2.fasta')


print(len(genome))



# What is the longest sequence and what is the shortest sequence? Is there more than one longest or shortest sequence? What are their identifiers?

mins = []
maxs = []
min = 100000
max = 0
for k,v in genome.items():
	l = len(v)
	if l < min:
		min = l
		mins = [k];
	elif l == min:
		mins+=k
		
	if l > max:
		max = l
		maxs = [k];
	elif l == max:
		maxs+=k
	
print ("MINS",mins)	
print ("MAXS",maxs)
print ("MIN",min)
print ("MAX",max)
#import collections
#od = collections.OrderedDict(sorted(genome.items()))

print(sorted(len(genome[i]) for i in genome.keys()))



'''
 In molecular biology, a reading frame is a way of dividing the DNA sequence of nucleotides into a set of consecutive, non-overlapping triplets (or codons). Depending on where we start, there are six possible reading frames: three in the forward (5' to 3') direction and three in the reverse (3' to 5'). For instance, the three possible forward reading frames for the sequence AGGTGACACCGCAAGCCTTATATTAGC are:

AGG TGA CAC CGC AAG CCT TAT ATT AGC

A GGT GAC ACC GCA AGC CTT ATA TTA GC

AG GTG ACA CCG CAA GCC TTA TAT TAG C

These are called reading frames 1, 2, and 3 respectively. An open reading frame (ORF) is the part of a reading frame that has the potential to encode a protein. It starts with a start codon (ATG), and ends with a stop codon (TAA, TAG or TGA). For instance, ATGAAATAG is an ORF of length 9.

Given an input reading frame on the forward strand (1, 2, or 3) your program should be able to identify all ORFs present in each sequence of the FASTA file, and answer the following questions: what is the length of the longest ORF in the file? What is the identifier of the sequence containing the longest ORF? For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier? What is the starting position of the longest ORF in the sequence that contains it? The position should indicate the character number in the sequence. For instance, the following ORF in reading frame 1:

>sequence1

ATGCCCTAG

starts at position 1.

Note that because the following sequence:

>sequence2

ATGAAAAAA

does not have any stop codon in reading frame 1, we do not consider it to be an ORF in reading frame 1.



!!!!
What is the length of the longest ORF appearing in reading frame 2 of any of the sequences?


'''


def min_positive(val1, val2):
	val = 0
	if val1 != -1 and val1 < val2:
		val = val1
	elif val2 != -1 and val2 < val1:
		val = val2
	return val

	
def get_orf_map(list, offset):
	start = 0
	map = {}
	for i in range(len(list)):
		if list[i] == 'ATG' and start == 0:
			start = i
		elif (list[i] == 'TAA' or list[i] == 'TAG' or list[i] == 'TGA') and start != 0:
			map[start*3 + offset+1] = (i - start + 1)*3
			start = 0
	return map	
	
def get_orf_maps(read):
	list1 = []
	list2 = []
	list3 = []

	for i in range(0, len(read)-2, 3):
		list1.append(read[i:i+3])

		list2.append(read[i+1:i+4])
		
		list3.append(read[i+2:i+5])
	
	if len(list2[-1]) < 3:
		del list2[-1]
	
	if len(list3[-1]) < 3:
		del list3[-1]
	
	#print (list1)
	#print (list2)
	#print (list3)
	return get_orf_map(list1, 0), get_orf_map(list2, 1), get_orf_map(list3, 2)
		
		
def deprecated_get_orf_map(read):
	orf_map = {}
	start = read.find('ATG')
	while start != -1:
		end = 100000
		end = min_positive(end, read.find('TAA', start + 2))
		end = min_positive(end, read.find('TAG', start + 2))
		end = min_positive(end, read.find('TGA', start + 2))
		if end > 0:
			end += 3
			orf_map[start] = end - start
			#print (read[start:end])
			start = read.find('ATG', end)
		else:
			break
	return orf_map
	
	
max = 0
id = ''
for k,v in genome.items():
	orf_map1, orf_map2, orf_map3 = get_orf_maps(v);
		
	if len(orf_map1) > 0: 
		read_max = sorted(orf_map1.values())[-1]
		if read_max > max:
			max = read_max
			id = k
		
	if len(orf_map2) > 0: 
		read_max = sorted(orf_map2.values())[-1]
		if read_max > max:
			max = read_max
			id = k
		
	if len(orf_map3) > 0: 
		read_max = sorted(orf_map3.values())[-1]
		if read_max > max:
			max = read_max
			id = k

print ("Longest ORF length", max)
print ("Longest ORF id", id)


max = 0
pos = -1
for k,v in genome.items():
	_, _ ,orf_map = get_orf_maps(v);
	for k1,v1 in orf_map.items():
		if v1 > max:	
			max = v1
			pos = k1
		
print("Longest ORF pos from 3", pos)
		
max = 0
orf_map1, orf_map2, orf_map3 = get_orf_maps(genome['gi|142022655|gb|EQ086233.1|16 marine metagenome JCVI_SCAF_1096627390048 genomic scaffold, whole genome shotgun sequence'])     #genome['gi|142022655|gb|EQ086233.1|101 marine metagenome JCVI_SCAF_1096627390048 genomic scaffold, whole genome shotgun sequence']);

if len(orf_map1)>0:
	lmax = sorted(orf_map1.values())[-1]
	if lmax > max:
		max = lmax

if len(orf_map2)>0: 
	lmax = sorted(orf_map2.values())[-1]
	if lmax > max:
		max = lmax

if len(orf_map2)>0:
	lmax = sorted(orf_map3.values())[-1]
	if lmax > max:
		max = lmax
				
print ("Longest ORF length in certain", max)



'''
A repeat is a substring of a DNA sequence that occurs in multiple copies (more than one) somewhere in the sequence. Although repeats can occur on both the forward and reverse strands of the DNA sequence, we will only consider repeats on the forward strand here. Also we will allow repeats to overlap themselves. For example, the sequence ACACA contains two copies of the sequence ACA - once at position 1 (index 0 in Python), and once at position 3. Given a length n, your program should be able to identify all repeats of length n in all sequences in the FASTA file. Your program should also determine how many times each repeat occurs in the file, and which is the most frequent repeat of a given length.
'''

frame_len = 7 #12
mega_dict = {}

max = 0
mega_frame = ''
mega_counter = 0 
#list = ['CAGTGCACAC', 'CAGACACAGACAGCAG']
for v in genome.values():
	for i in range(len(v) - frame_len + 1):
		frame = v[i:i+frame_len]
		#print(frame)
		if frame in mega_dict:
			mega_dict[frame] += 1
			if mega_dict[frame] > max:
				max = mega_dict[frame]
				mega_frame = frame
				mega_counter = 0
				#print(frame);
			elif mega_dict[frame] == max:
				mega_counter+=1
		else:
			mega_dict[frame] = 1


print("MAX Repeat", mega_frame)
print("MAX Repeat num", max)
print("Number of max repeats of the same freq", mega_counter+1)
print(mega_dict['GCGGCCG'])			
print(mega_dict['CGCGCCG'])			
print(mega_dict['TGCGCGC'])			
print(mega_dict['CATCGCC'])			