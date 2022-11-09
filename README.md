# README
Major components include:
* graph_decompose: subgraph decomposing logic
* graph_counting: sampling to count logic
* peregrine: a previous enumerating work as baseline

# Graph Counting
Count number of certain patterns in a large graph

## Quick Test
```
cd graph_counting/src/{branch}
make clean ; make all
mpirun -n 4 -f host ./GraphCounting.out graphs/mico/mico.undigraph patterns/triangle 40000000 40 4
```

Output:
```
raw results after merge nodes (total_prob_inverse, total_sampled_times): 502464449642959.8 40000000
pattern_nums time_consumed(us) sampling_times sampling/second
12561611.24107399 2065807 40000000 19362893.04857617
```

## Usage:
```
GraphCounting graph_file subgraph_file sampling_times(separated by \',\' 10000,20000) #threads(separated by \',\' 1,4) #processes
```
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


## Memcached - Version (batched_distributed_samplers and fail_prob_aware_distributed_samplers branch)
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



## MPI Installation
We use MPI and multi-threading for load balancer of the replicated graph setting.

First set up ssh between node0 and other nodes.

MPI installalation:

Download MPICH3 (http://www.mpich.org/downloads/) from mpich.org and extract the contents of the MPICH package to some temporary location, before compiling the source code to binaries. We download MPICH v3.2 (http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz).
```
cd {the-directory-containing-downloaded-mpich-package}
tar -xvzf mpich-3.2.1.tar.gz
cd mpich-3.2.1
```
Choose a directory for installing MPICH3 (we use /usr/local/mpich), and then compile and install MPICH3.
```
./configure -prefix=/usr/local/mpich --disable-fortran
sudo make
sudo make install
```
Append the following two environment variables to the file $HOME/.bashrc. Here, $HOME, ~ and /home/{your_username} are equivalent, which is your home folder.
```
export MPICH_HOME=/usr/local/mpich
export PATH=$PATH:$MPICH_HOME/bin
```
Compile the file with the command source **$HOME/.bashrc**.
