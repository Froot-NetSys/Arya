# input num_nodes and num_edges of the graph you wish to test

python ../generate_random_graph.py $1 $2 test ; python ../graph_decomp.py test ; cat test

# Use webgraphviz to test: graph{1--2 .....}
