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
# input : sequence from datafile from commandline
# output : probabilistic seq/matrix
# this is also to show Colby the structure of a dictionary
# usually I just use standard dictionary but it is saved in non-ordered way. Just to show you, using an ordered dictionary
def parse_seq(data):
        prob_seq = collections.OrderedDict()
        with open(data) as f:
                for seq in f:
                        seq = seq.strip()
                        for nuc in nucleotides:
                                for i in range(0,len(seq)):
                                        if nuc == seq[i]:
                                                prob_seq[nuc,i] = 1
                                        else:
                                                prob_seq[nuc,i] = 0
                return prob_seq

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
                                print "File not found !"
                                help()
                                sys.exit(1)
        try:                    
                parse_seq(file_name)
        except NameError:
                print "Need to specify a file !"
                help()
                sys.exit(1)
