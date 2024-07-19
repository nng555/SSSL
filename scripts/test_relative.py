from kmer_trie import KmerTrie

test = 'ACTGATCGTAGCATGTGGGGCTACTCCGCGGCGCGCATTACGACCG'
etrie = KmerTrie(8)
etrie.insert(test, 'r1')
etrie.insert(test, 'r2')
etrie.relative_entropy('r1', 'r2', 8)
