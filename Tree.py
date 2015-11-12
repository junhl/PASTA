class Tree(object):

    def __init__(self, root):
        """
        root : the root of the tree
        nodes: a list of all the nodes in the tree
        """
        self.root = root
        

    def harvest(self):
        """ Returns a list of strings, all selected k mers in the subtree specified by root, contained in the leaves of
        the tree
        """
        
        #Perform a depth first search on the graph. If the currently selected node is a leaf, add it to a list of them
        leaves = [];
        stack = []
        stack.append(self.root)
        max_depth = len(self.root.part_k_mer)
        #While the stack is not empty
        while stack:
            current_node = stack.pop()
            for child in current_node.children:
                #If the child is a leaf
                if '-' not in child.part_k_mer:
                    leaves.append(child.part_k_mer)
                else:
                    stack.append(child)
        return leaves

