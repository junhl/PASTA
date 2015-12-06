# Authors  Jun Ho Lee, Colby Simpson
# PASTA - Probabilistic Alignment Search Tool Advance
# Python Tree Implementation

class Node(object):

    def __init__(self, parent, children, prob_score, prob_dist, part_k_mer, depth):
        """
        children: a list of child nodes
        prob_score: the currently calculated probability of our probabilistic string
                    emitting the partial k mer as a substring
        prob_dist. The probability distribution of the ith letter of our probabilistic
                   string. All nodes at depth i from the root will have the same distribution
        part_k_mer: A string representing the partially constructed k mer
        depth: an integer representing a nodes depth from the root.
        """
        self.parent = parent
        self.children = children
        self.prob_score = prob_score
        self.prob_dist = prob_dist
        self.part_k_mer = part_k_mer
        self.depth = depth
      
   
