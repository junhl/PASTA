# Authors  Jun Ho Lee, Colby Simpson
# PASTA - Probabilistic Alignment Search Tool Advance
# Database Builder
'''
   The database builder parses the given nucleotide and probability sequences
   into a workable data structrue which is then saved for use in later calculations
'''


probSequence ={}
for i in range(len(n)):
        """ Take the probability, match it to the given nucleotide in
            the dictionary at sequence index i. Then calcualte the
            probability of other nucleotides
        """
        givenProb = p[i]
        for nuc in ['A', 'C', 'G', 'T']:
                if nuc == n[i]:
                        probSequence[(i, nuc)] = float(p[i])
                else:
                        probSequence[(i, nuc)] = (1 - float(givenProb))/3.0

pickle.dump(probSequence, open("database.dat", "w+"))
