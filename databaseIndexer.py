# Authors  Jun Ho Lee, Colby Simpson
# PASTA - Probabilistic Alignment Search Tool Advance
# Database Indexer

"""
Given a probabilistic sequence S, we wish to find all deterministic k - mers of S
emitted with probability greater than a threshold value t. We return a dictionary of the form
dict[index] = [list of all k mers with instantiation probability greater than t starting at this index]
"""
from Node import Node
from Tree import Tree
from math import log

def database_indexer(k,threshold, nucleotides, probabilities):
    """
    k : length of k - mer
    t : threshold value
    nucleotides : string which is the filename of our list of nucleotides
    prob : string which is the filename of our list of probabilites
    """
    t = log(threshold)
    nuc = []
    #The probabilities in string format.
    prob_string = []
    probable_emissions = {}
    
    
    with open (nucleotides) as f:
        for line in f:
            for char in line:
                nuc.append(char)
    
    with open(probabilities) as f:
        for line in f:
            prob_string = line.split(' ')
    #Convert the probability values from strings to floating point numbers
    prob = [float(i) for i in prob_string]
    
    distributions = []
    for i in range(len(nuc)):
        distributions.append(get_ith_distribution(i, nuc, prob))
    #A list which has as the element at index i the log of the maximum probability of an instantiated suffix of S starting at index i. We currently don't make any calls to this
    max_word_probability = []
    for i in range(len(distributions)):
        max_word_probability.insert(i, prob_upper_bound(distributions, i, k))
                      
    for start_idx in range(len(nuc)- k +1):
        #This is "window" of the probabilistic sequence of length k from which we will look for k - mers instantiating with probability > t
        distributions_window = distributions[start_idx:start_idx + k ]
        """
        Set up the root of the k mer tree. It currently has no children, a probability of 1, the distribution
        of the first letter of the probabilistic sequence S, and a string of k"-" characters, showing that
        all letters in the k - mer are not instantiated
        """
        base_string = list('-'*k)
    
        # Look for characters in S which are certain. This will help us minimize the number of nodes created
        for i in range(len(distributions_window)):
            if distributions[i][0][1] == 1:
                base_string[i] = distributions_window[i][0][0]

        root_distro = distributions_window[0]
        root = Node(None, [], 0, root_distro, base_string, 0)
        #Set up the tree using Gardener
        gardener(t, root, 0, k, distributions_window)
        #Once nodes have been linked, consolidate to a tree
        k_mer_tree = Tree(root)
        probable_emissions[start_idx] = k_mer_tree.harvest()
    
    return probable_emissions


def gardener(t,seed, current_depth, max_depth, distributions):
    """ --GROWS THE TREE LOL--
    Given a node seed, which represents a partially generated k - mer of length cur_Idx, set up
    its children to represent all partially generated k - mers of length cur_Idx + 1
    Then recurse onto the child nodes, returning in the end a tree whose leaves are those k mers which are instantiated with probability > threshold. 
    """
    
    """
    rule of growth: every path in the tree starting from the root is of length <= (length of k gram).
     -    those of length < k are branches which terminate because they fall below the acceptable probability threshold
     -    those of length = k have all characters in their k -grams instantiated
    """
    if seed.depth < max_depth:
        # If the current letter of the k - mer is not instantiated/
        if seed.part_k_mer[current_depth] == '-' :
            #pair  = ('N, p)
            for pair in seed.prob_dist:
                # some test variables to make debugging easier
                # current_probability = seed.prob_score + log(pair[1])
                # upper_bound = prob_upper_bound
                # test = current_probability + upper_bound
                             
                if seed.prob_score + log(pair[1])  > t:
                    #Using list() ensures that we get a copy of the list, and won't modify the k-mer stored in the Seed node
                    new_k_mer = list(seed.part_k_mer)
                    new_k_mer[current_depth] = pair[0]
                    #Leaves don't generate children, so they don't get a distribution
                    if seed.depth == max_depth - 1:
                        child = Node(seed,[], seed.prob_score + log(pair[1]), None, new_k_mer, current_depth + 1)
                    else:
                        child = Node(seed,[], seed.prob_score + log(pair[1]), distributions[current_depth +1], new_k_mer, current_depth+1)
                    seed.children.append(child)
        # else, the current node must have a certain nucleotide at the current position
        else:
            #If the last letter in the k - mer is certain, then we can stop growth here
            if seed.depth == max_depth - 1:
                pass
            else:
                """Check the seed probability score here?"""
                child = Node(seed,[], seed.prob_score, distributions[current_depth + 1], seed.part_k_mer, current_depth + 1)
                seed.children.append(child)                           
        # Now grow each of the child nodes
        for child in seed.children:
            gardener(t,child, current_depth + 1, max_depth, distributions) 
        

def get_ith_distribution(i, nuc, prob):
    """
    Return a list representing the probability distribution of the ith letter of the
    probabilistic sequence. The list is sorted in order of largest probability to least
    """
    distro = []
    distro.append( (nuc[i], float(prob[i])) )
    if prob[i] == 1.0:
        return distro
    else:
        available_nucs = [x for x in ['A', 'C', 'T', 'G'] if x != nuc[i]]
        for char in available_nucs:
            distro.append( (char, (1-float(prob[i]))/3.0))
        """This sorting may not be neccessary by the construction of the list given, but might be a good idea to keep it in anyways since it affords some robustness to the input we can take"""
        distro.sort(key=lambda tup: 1/tup[1])
        # Is it worth it? Let me work it. I put my list down, flip it and reverse it
        #distro_sorted = distro[::-1]
    #return distro_sorted
    return distro

"""
As of now we don't make any calls to this method. However at a later time we might be able to make use of it for a more effective branch and bound strategy

"""
def prob_upper_bound(distributions, i, k):
    """ calculate the log of the maximum instantiation probability of all suffixes of S starting at i with length k"""
    ub = 0;
    window = distributions[i:i+k]
    for j in range(len(window)):
        # distributions[i][0][1] is the maximum probability of a nucleotide appearing at position i of S
        ub += log(window[j][0][1])
    return ub
