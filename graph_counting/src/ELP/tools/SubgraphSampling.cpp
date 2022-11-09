#include <iostream>
#include <string.h>
#include <chrono>
#include <atomic>
#include <thread>
#include <pthread.h>
#include <mutex>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iomanip>
#include "graph.h"
#include <ctime>

using namespace std;

 // Subgraph Sampling

    // uniform_edge_sampling: sampling edges with probablity sampling_ratio / 100.
    void uniform_edge_sampling(Graph & large_graph, string subgraph_file_, uint32_t sampling_ratio) // sampling_ratio%
    {
        srand(time(NULL));
        std::ofstream outfile(subgraph_file_);
        uint32_t edge_count = 0;
        for (uint32_t i = 0; i < large_graph.total_edges_; i++)
        {
            if (large_graph.edge_list_[i].v_start > large_graph.edge_list_[i].v_end) continue;
            if (std::rand() % 100 + 1 > sampling_ratio) continue;
            outfile << large_graph.edge_list_[i].v_start << " " << large_graph.edge_list_[i].v_end << endl;
            edge_count++;
        }
        cout << "total_edges = " << large_graph.total_edges_ << endl;
        cout << "edge_count = " << edge_count << endl;
    }
    

 // uniform_node_sampling: sampling nodes with probablity sampling_ratio / 100, and obtain all neighbors of each node.
    void uniform_node_sampling(Graph & large_graph, string subgraph_file_, uint32_t sampling_ratio) // sampling_ratio%
    {   
        std::srand(static_cast<unsigned int>(std::time(nullptr))); 
        std::ofstream outfile(subgraph_file_);
        uint32_t node_count = 0;
        for (auto it = large_graph.v_table_.begin(); it != large_graph.v_table_.end(); it++)
        {
            if (std::rand() % 100 + 1 > sampling_ratio) continue;
            node_count++;
            for (uint32_t idx = 0; idx < it->second.degree; idx++) {
                outfile << it->first << " " << large_graph.edge_list_[it->second.edge_index + idx].v_end << endl;
            }
        }
        cout << "total_nodes = " << large_graph.v_table_.size() << endl;
        cout << "node_count = " << node_count << endl;
    }

int main(int argc, char *argv[])
{
    // Parse running parameters
    if (argc < 4) {
        cout << "Sampling subgraphs in a large graph \n"
             << "Usage: SubgraphSampling graph_file_path edge_subgraph_file_path sampling_ratio(%)"
             << endl;
        return 1;
    }
    string graph_file_path = argv[1];
    string subgraph_file_path_edge = argv[2];
    // string subgraph_file_path_node = argv[3];
    uint32_t sampling_ratio = stoi(argv[3]); // sampling_ratio% = sampling_ratio / 100
    Graph large_graph(graph_file_path);
    large_graph.load();

    uniform_edge_sampling(large_graph, subgraph_file_path_edge, sampling_ratio);
    cout << "* uniform edge sampling finished!" << endl;
    // uniform_node_sampling(large_graph, subgraph_file_path_node, sampling_ratio);
    // cout << "* uniform node sampling finished!" << endl;
}