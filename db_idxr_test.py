# database indexer debug script
from databaseIndexer import database_indexer

n = 'indexertext_nuc.txt'
p = 'indexertext_prob.txt'
k = 5
d = database_indexer(3, 0.054,n,p)
print d
