# Graph Counting
Count number of certain patterns in a large graph

## Quick Test
```
make clean ; make
./GraphCounting.out ./test_graphs/test_10_nodes ./patterns/5_cycle 10000

output:
.............
Estimating finds 6.0928 patterns
pattern_nums time_consumed(us)
6.0928 52056
```

The example counts number of 5-cycles in test_10_nodes small figure, with 10000 sampling times.

## Usage:
GraphCounting.out graph_file subgraph_file sampling_times

* graph_file: edge list of the large graph ; sorted by vertex ids, each edge has two copies (e.g., (1,2), (2,1)), example: 
```
1 2
1 3
2 1
3 1
```
* subgraph_file: pattern to count ; example: ./patterns/triangle_2_star
    * the subgraph_file is supposed to be output from the graph_decomposing logic
```
# first line: num of odd-cycles (k), num of stars (j) in the pattern
# following k lines: odd-cycles
# following j lines: stars (the first vertex is the center)
# following lines: remaining edges to test
1 1   # indicates 1 odd-cycle, 1 star
1 2 3 # the first and only odd-cycle
4 5 6 # the first and only star
2 4   # only one remaining edge
```


## Memcached - Version (distributed_estimators branch)
* Step 0: Install Memcached and libmemcached
    * `cd scripts ; bash config.sh`
    * remote config: `cd scripts; bash remote_config.sh node1` to config node1
* (Not Used) Start Memcached with: 
    * enable remote connect: `pkill memcached ; ufw allow 11211; memcached -d -m 100000 -p 11211 -u root -M`
    * NOTE: no need to start with automatic loading
* Step 1: load graph:
    * `cd scripts ; python multi_machine_load.py graph_directory`
        * e.g. `python multi_machine_load.py ../split_graph/mico`
        * single machine load: `python single_machine_load.py graph_path edge_start_id`, e.g. `python single_machine_load.py global_subgraph_0 0`
    * in the graph_directory:
        * machine_info: must contain num_machines and total_edges in the first line
        * each line is a memcached node
* Step 2: distributed estimators:
    * `python distributed_estimators.py machine_info pattern total_num_estimators num_threads_per_node`
    * e.g. `python ./distributed_estimators.py ./split_graph/mico/machine_info ./patterns/5_cycle_2_star 10000000 20`

