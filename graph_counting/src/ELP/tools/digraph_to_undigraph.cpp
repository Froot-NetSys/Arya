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
#include <unordered_map>
#include <set>
#include <fstream>

using namespace std;

unordered_map<uint32_t, set<uint32_t> > edge_list;

int main(int argc, char *argv[])
{
    // Parse running parameters
    if (argc < 3) {
        cout << "generate undigraph from digraph \n"
             << "Usage: Undigraph.out digraph_file_path undigraph_file_path"
             << endl;
        return 1;
    }

    string digraph_file_path = argv[1];;
    string undigraph_file_path = argv[2];
    std::ifstream infile(digraph_file_path);
    std::ofstream outfile(undigraph_file_path);
    
    uint32_t v_start, v_end;
    while (infile >> v_start >> v_end)
    {   
        if (edge_list.find(v_start) == edge_list.end())
        {
            set<uint32_t> adj;
            adj.insert(v_end);
            edge_list.emplace(v_start, adj);
        }
        else
        {
            edge_list.at(v_start).insert(v_end);
        }
        if (edge_list.find(v_end) == edge_list.end())
        {
            set<uint32_t> adj;
            adj.insert(v_start);
            edge_list.emplace(v_end, adj);
        }
        else
        {
            edge_list.at(v_end).insert(v_start);
        }
    }
    cout << "finish loading digraph" << endl;
    for (auto it = edge_list.begin(); it != edge_list.end(); it++)
    {
        for (auto adj = it->second.begin(); adj != it->second.end(); adj++)
            outfile << it->first << " " << *adj << endl;
    }
    cout << "finish writing undigraph" << endl;
}