import sys
import os

def load_and_split_graph(graph_in, result_graph_dir, partition_number):
    # load in graph
    edges = []
    nodes = {}
    degree = {}
    for line in graph_in:
        if '#' in line:
            continue
        line_split = line.strip('\n').split(' ')
        start_vertex = int(line_split[0])
        end_vertex = int(line_split[1])
        
        nodes[start_vertex] = 1
        nodes[end_vertex] = 1
        
        edges.append([start_vertex, end_vertex])
        if start_vertex not in degree:
            degree[start_vertex] = 0
        degree[start_vertex] += 1

    nodes_list = nodes.keys()
    node_to_index = {}
    for i in range(len(nodes_list)):
        node_to_index[nodes_list[i]] = i

    n = len(nodes_list)
    m = len(edges)
    print n, m

    sum_degree = {}
    last = 0
    sum_degree[0] = degree[nodes_list[0]]
    for i in range(1, len(nodes_list)):
        sum_degree[i] = degree[nodes_list[i]] + sum_degree[i-1]

    # debug
    #for i in range(n):
    #    print nodes_list[i], degree[nodes_list[i]], sum_degree[i]


    # split graph into subgraphs with 2 versions
    # assume edges are sorted by vertex

    # consider discrete node indexes to correctly calculate indexes of subgraphs
    subgraph_node_index = []
    subgraph_edge_index = []
    begin_node_index = 0
    begin_edge_index = 0
    for i in range(n % partition_number):
        subgraph_node_index.append([begin_node_index, begin_node_index + n / partition_number])
        if begin_node_index == 0:
            delta_edge_index = sum_degree[n / partition_number]
        else:
            delta_edge_index = sum_degree[begin_node_index + n / partition_number] - sum_degree[begin_node_index - 1]
        subgraph_edge_index.append([begin_edge_index, begin_edge_index + delta_edge_index - 1])

        begin_node_index += n / partition_number + 1
        begin_edge_index += delta_edge_index

    for i in range(n % partition_number, partition_number):
        subgraph_node_index.append([begin_node_index, begin_node_index + n / partition_number - 1])
        if begin_node_index == 0:
            delta_edge_index = sum_degree[n / partition_number - 1]
        else:
            delta_edge_index = sum_degree[begin_node_index + n / partition_number - 1] - sum_degree[begin_node_index - 1]
        subgraph_edge_index.append([begin_edge_index, begin_edge_index + delta_edge_index - 1])

        begin_node_index += n / partition_number
        begin_edge_index += delta_edge_index

    # global sampling version: create edge files of subgraphs
    with open(result_graph_dir + "machine_info", "w") as info_out:
        info_out.write(str(partition_number) + "\n")
        for i in range(partition_number):
            info_out.write(node_ip[i] + " ")
            info_out.write(str(nodes_list[subgraph_node_index[i][0]]) + " " + str(nodes_list[subgraph_node_index[i][1]]) + " ")
            info_out.write(str(subgraph_edge_index[i][0]) + " " + str(subgraph_edge_index[i][1]) + '\n')
            with open(result_graph_dir + "global_subgraph_" + str(i), "w") as graph_out:
                for j in range(subgraph_edge_index[i][0], subgraph_edge_index[i][1] + 1):
                    graph_out.write(str(edges[j][0]) + " " + str(edges[j][1]) + "\n")
    
    # local sampling version: create edge files of subgraphs
    for i in range(partition_number):
        with open(result_graph_dir + "local_subgraph_" + str(i), "w") as graph_out:
            # graph_out.write("# nodes: " + str(subgraph_node_index[i][0]) + "-" + str(subgraph_node_index[i][1]) + "\n")
            # graph_out.write("# edges: " + str(subgraph_edge_index[i][0]) + "-" + str(subgraph_edge_index[i][1]) + "\n")
            print "local_sampling_subgraph " + str(i)
            # print "# nodes: " + str(subgraph_node_index[i][0]) + "-" + str(subgraph_node_index[i][1])
            # print "# edges: " + str(subgraph_edge_index[i][0]) + "-" + str(subgraph_edge_index[i][1]) + '\n'
            for j in range(subgraph_edge_index[i][0], subgraph_edge_index[i][1] + 1):
                if node_to_index[edges[j][1]] >= subgraph_node_index[i][0] and node_to_index[edges[j][1]] <= subgraph_node_index[i][1]:
                    graph_out.write(str(edges[j][0]) + " " + str(edges[j][1]) + "\n")
                
        
if __name__ == "__main__":
    print "python split_graph.py original_graph_file split_graph_dir partition_number"
    original_graph_file = sys.argv[1]
    split_graph_dir = sys.argv[2]
    partition_number = int(sys.argv[3])
    node_ip = ["10.10.1.1", "10.10.1.2", "10.10.1.3", "10.10.1.4", "10.10.1.5", "10.10.1.6", "10.10.1.7", "10.10.1.8", "10.10.1.9", "10.10.1.10", "10.10.1.11", "10.10.1.12", "10.10.1.13", "10.10.1.14", "10.10.1.15", "10.10.1.16"]
    with open(original_graph_file, "r") as graph_in:
        load_and_split_graph(graph_in, split_graph_dir, partition_number)
