import sys
import random



if len(sys.argv) < 4:
    print("To use: python generate_random_graph.py num_nodes num_edges graph_output_file_name")
    print("we will generate random unique undirected edges with nodes from [0, num_nodes)")
    exit(1)

num_nodes = int(sys.argv[1])
num_edges = int(sys.argv[2])
graph_file = sys.argv[3]

edges = set()

with open(graph_file, 'w') as graph_output:
    for i in range(num_edges):
        src_node = -1
        target_node = -1
        while True:
            src_node = random.randint(0, num_nodes-1) 
            target_node = random.randint(0, num_nodes-1)
            if src_node == target_node:
                continue
            if src_node < target_node:
                if str(src_node) + '->' + str(target_node) not in edges:
                    break
            else:
                if str(target_node) + '->' + str(src_node) not in edges:
                    break

        edges.add((src_node, target_node))
        edges.add((target_node, src_node))
        
    for edge in sorted(edges, key=lambda tup: (tup[0], tup[1])):
        graph_output.write(str(edge[0]) + ' ' + str(edge[1]) + '\n')   

