# Authors  Jun Ho Lee, Colby Simpson
# PASTA - Probabilistic Alignment Search Tool Advance

import collections
import pickle

n = []
p = []
with open ('chr22.fa') as f:
	for line in f:
		for nuc in line:
			n.append(nuc)
print len(n)

with open ('chr22.conf') as f:
	for line in f:
		p = line.split(' ')
print len(p)
for x in range(10):
	print p[x]


probSequence ={}
for i in range(len(n)):
        """
        Take the probability, match it to the given nucleotide in
        the dictionary at sequence index i. Then calculate the
        probability of other nucleotides
        """
        givenProb = p[i]
        for nuc in ['A', 'C', 'G', 'T']:
                if nuc == n[i]:
                        probSequence[(i, nuc)] = float(p[i])
                else:
                        probSequence[(i, nuc)] = (1 - float(givenProb))/3.0

pickle.dump(probSequence, open("database.dat", "w+"))
