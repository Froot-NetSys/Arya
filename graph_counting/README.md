# Graph Counting
Count number of certain patterns in a large graph

## Quick Test
```
cd src/{branch}
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
