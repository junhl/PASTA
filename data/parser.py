# Authors  Jun Ho Lee, Colby Simpson
# PASTA - Probabilistic Alignment Search Tool Advance

import collections

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