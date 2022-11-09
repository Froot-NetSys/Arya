# Graph Decomposition
Decompose a pattern graph into odd cycles and stars

# Installation:
1. Python 2.7
2. Install scipy: 

```
sudo apt-get install build-essential gfortran libatlas-base-dev python-pip python-dev
sudo pip install --upgrade pip
sudo pip install numpy
sudo pip install scipy
```

# Main program:
graph_decomp.py

## Quick Test
```
cd ./test_graph
./test.sh 10 20    // 10 is #nodes, 20 is #edges
```

Will random generate a graph pattern with assigned #nodes and #edges, and then decompose it.

## Usage:
Input: undirected graph(edges), example: test_graph/test
Usage: `python graph_decomp.py test`
Output: see print out 
Example:
```
The optimal fractional cover:  [0.5 0.5 0.5 0.5 0.5 0.  1.  1.  0.  0.  1. ]
0 -- 1 [label= 0.5 ,color= green ]
1 -- 2 [label= 0.5 ,color= green ]
2 -- 3 [label= 0.5 ,color= green ]
3 -- 4 [label= 0.5 ,color= green ]
0 -- 4 [label= 0.5 ,color= green ]
4 -- 7 
5 -- 7 [label= 1.0 ,color= red ]
6 -- 7 [label= 1.0 ,color= red ]
7 -- 8 
7 -- 9 
8 -- 9 [label= 1.0 ,color= red ]
==== decomposing the support of the fractional edge cover ====
Following are odd, vertex-disjoint cycles within the graph
[[0, 4, 3, 2, 1, 0]]
Following are stars(connected components)
[[8, 9], [5, 7, 6]]

```

The first part is the graph with labels for values and colors to identify it belongs to a cycle or star.

The scond part is the cycles and stars.

## Visualization:
Use: http://www.webgraphviz.com
Example graph code:

```
graph G {
  6 -- 9 [label= 0.5 ,color= green ]
  3 -- 8 [label= 1.0 ,color= red ]
  1 -- 2 [label= 0.5 ,color= green ]
  1 -- 5 [label= 0.5 ,color= green ]
  5 -- 9 [label= 0.5 ,color= green ]
  7 -- 8 
  4 -- 7 [label= 1.0 ,color= red ]
  2 -- 6 [label= 0.5 ,color= green ]
  2 -- 3 
  2 -- 7 
}

```
