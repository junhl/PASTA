# Authors  Jun Ho Lee, Colby Simpson
# PASTA - Probabilistic Alignment Search Tool Advance

import collections
import sys
import os
import pickle

db_path = os.path.relpath("data/database.dat")
nucleotides = ['A','C','G','T']

db = pickle.load(open (db_path, 'r'))
# a helper function that may be useful in later on. 
def parse_seq(data):
    prob_seq = [] #collections.OrderedDict()
    with open(data) as f:
            for seq in f:
				seq = seq.strip()
				for i in range(0,len(seq)):
					for nuc in nucleotides:
						if nuc == seq[i]:
							prob_seq.append(nuc)
						'''else:
							prob_seq[i,nuc] = 0.00'''
    return prob_seq
'''
print db[(0,'C')]
print db[(0,'T')]

print db[(68,'A')]
print db[(69,'A')]
print db[(70,'A')]


print db[(156,'A')]
print db[(157,'A')]
print db[(158,'A')]

print db[(604316,'A')]
print db[(604317,'A')]
print db[(604318,'A')]
'''
def makeWords(db_dic, query):
	candidates = []
	positions = []
	for j in range(0,len(query)-2):
		q = query[j]+query[j+1]+query[j+2]
		candidates.append(q)
		positions.append((j,j+2))
	words = []
	for i in range(0,-2+(len(db_dic)/4)):
		for n1 in nucleotides:		
			if db_dic[(i,n1)] == 1.00:
				for n2 in nucleotides:
					if db_dic[(i+1,n2)] == 1.00:
						for n3 in nucleotides:
							if db_dic[(i+2,n3)] == 1.00:
								w = n1+n2+n3
								if w in candidates:
									indices = [k for k, x in enumerate(candidates) if x == w]
									for x in indices:
										words.append((w,(i,i+2),positions[x]))
	'''print len(words)
	print words[0]
	print words[1]
	print words[2]
	print words[-2]
	print words[-1]'''
	# words is a list of (seed/word, start and end of word in db, start and end of word in query)
	return words

def extend(seed, query, threshold): # ignore the treshold for now
	scores = []
	scores2 = []
	for i in range (len(seed)): # extend for each case
		numerator = 3.0 # see notice.txt for the reasoning behind this
		denominator = 3.0
		
		first = seed[i][2][0]
		end = seed[i][2][1]
		#print first,end, len(query), (first != 0), (end != len(query)-1)
		while first != 0 or end != len(query)-1: # extend until both ends are reached to make it false false
			if first != 0: # if not reached end of left, extend to the left
				nuc = query[first-1] # the nucleotide that we are adding to alignment
				denominator += 1.0
				first -= 1
				
				
				if (db[seed[i][1][1]+1, nuc] >= threshold):
					numerator += db[seed[i][1][0]-1, nuc]
					seed[i] = (nuc + seed[i][0], (seed[i][1][0]-1, seed[i][1][1]), (seed[i][2][0]-1, seed[i][2][1]))
				
				else:
					numerator += threshold
					seed[i] = (nuc + seed[i][0], (seed[i][1][0]-1, seed[i][1][1]), (seed[i][2][0]-1, seed[i][2][1]))
				
			if end != len(query)-1: # if not reach the end of sequence, extend to the right. Notice that its not elif
				nuc = query[end+1] # the nucleotide that we are adding to alignment
				denominator += 1.0
				end += 1
				
				if (db[seed[i][1][1]+1, nuc] >= threshold):
					numerator += db[seed[i][1][1]+1, nuc]
					seed[i] = (seed[i][0]+nuc, (seed[i][1][0], seed[i][1][1]+1), (seed[i][2][0], seed[i][2][1]+1))
				
				else:
					numerator += threshold
					seed[i] = (seed[i][0]+'-', (seed[i][1][0], seed[i][1][1]+1), (seed[i][2][0], seed[i][2][1]+1))
							
		scores.append(numerator/denominator)			
		scores2.append(numerator)
	print scores.index(max(scores)), seed[scores.index(max(scores))]
	print 'Your query sequence', seed[scores.index(max(scores))][0], ' is best aligned with score of', scores[scores.index(max(scores))], 'at ', seed[scores.index(max(scores))][1]
	print 'Using method 2 : ', seed[scores2.index(max(scores2))][0], ' is best aligned with score of ', scores2[scores2.index(max(scores2))], 'at', seed[scores2.index(max(scores2))][1]
	return scores
# help for general users. To be expanded as we add more parameters
def help():
  print """
  Required:
    -d <file_path> 
      A file containing the query sequence

  e.g.
    python pasta.py -d query.txt
    """

 
if __name__ == "__main__":
	#Just a little hack for those working in IDLE and not off the command line
    #sys.argv = ["pasta.py", "-d query.txt"]
    param = sys.argv
    l = len(param)
    i = 1
       
    while i < l:
        options = param[i]
        if not options.startswith('-'):
                i += 1
        else:
            if options == '-d':  
                file_name = param[i+1]
                i += 2
            if not os.path.isfile(file_name):
                print 'File not found !'
                help()
                sys.exit(1)
	seq = parse_seq(file_name)
	print "query is ", seq
	print db[16216,'A']
	print db[16217,'A']
	print db[16218,'A']
	print db[16219,'G']
	print db[16220,'C']
	print db[16221,'C']
	print db[16622,'C']
	words = makeWords(db,seq)
	print len(words)
	
	extend(words, seq, 0.25)