import sys
from cvxopt import matrix, solvers
from scipy.optimize import linprog

comment_flag = False  # whether print detailed steps

class Node:
   def __init__(self, node_id):
       self.node_id = node_id
       self.neighbors = set()   # set of neighbor id

   def add_neighbor(self, neighbor):
       self.neighbors.add(neighbor)
       #self.neighbors.append(neighbor)

   def get_node_id(self):
       return self.node_id

   def get_neighbors(self):
       return self.neighbors

   def dump_node(self):
       print self.node_id, ':', self.neighbors

class Graph:
    def __init__(self): 
        self.nodes = {}         # node id to Node instance
        self.num_edges = 0
        self.edge_id_map = {}   # (src node id, target node id) to edge id
        self.edge_backmap = {}  # edge_id to (src node id, target node id)
    def add_edge(self, src_node_id, target_node_id):
        if comment_flag:
            print 'adding: ', src_node_id, target_node_id
        if src_node_id == target_node_id:
            print '=========== ALERT: same src and target nodes'
            exit(1)
        if src_node_id not in self.nodes:
            self.nodes[src_node_id] = Node(src_node_id) 
        if target_node_id not in self.nodes:
            self.nodes[target_node_id] = Node(target_node_id)

        self.nodes[src_node_id].add_neighbor(target_node_id)
        self.nodes[target_node_id].add_neighbor(src_node_id)

        # set id(the current num_edges) for this edge 
        if src_node_id < target_node_id:    # the samller one always comes first 
            edge_string = str(src_node_id) + '->' + str(target_node_id)
        else:
            edge_string = str(target_node_id) + '->' + str(src_node_id)
        if edge_string not in self.edge_id_map:
            self.edge_id_map[edge_string] = self.num_edges
            self.edge_backmap[self.num_edges] = edge_string
            self.num_edges += 1

    def get_edge_id(self, src_node_id, target_node_id):
        if src_node_id < target_node_id:    # the samller one always comes first 
            edge_string = str(src_node_id) + '->' + str(target_node_id)
        else:
            edge_string = str(target_node_id) + '->' + str(src_node_id)
        if edge_string not in self.edge_id_map:
            print 'ALERT: edge_string cannot find'
            exit(1)
        return self.edge_id_map[edge_string]
    
    def translate_edge_id(self, edge_id):
        return int(self.edge_backmap[edge_id].split('->')[0]), int(self.edge_backmap[edge_id].split('->')[1]) # the samller one always comes first 
    
    def get_nodes(self):
        return self.nodes.values() 
    
    def get_edges(self):
        return self.edge_id_map.values() 
    
    def count_nodes(self):
        return len(self.nodes)
    
    def count_edges(self):
        return self.num_edges

    def get_neighbors(self, node_id):
        return self.nodes[node_id].get_neighbors()

    # Method to retrieve connected components 
    # in an undirected graph
    def DFSUtil(self, temp, v, visited): 
  
        # Mark the current vertex as visited 
        visited[v] = True
  
        # Store the vertex to list 
        temp.append(v) 
  
        # Repeat for all vertices adjacent 
        # to this vertex v 
        for i in self.nodes[v].get_neighbors():
            if visited[i] == False: 
                  
                # Update the list 
                temp = self.DFSUtil(temp, i, visited) 
        return temp 
 
    def connectedComponents(self): 
        visited = {}
        cc = [] 
        for node in self.nodes.values():
            visited[node.get_node_id()] = False
        for node in self.nodes.values():
            v= node.get_node_id()
            if visited[v] == False:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
        return cc 

    def dump_graph(self):
        print '===== dump the graph: '
        for node in self.nodes:
            self.nodes[node].dump_node()
 


def parse_graph(graph_in):
    G = Graph()

    for line in graph_in:
        G.add_edge(int(line.split('--')[0]), int(line.split('--')[1]))    # this is the format of graphviz
        if comment_flag:
            G.dump_graph()
    return G
'''
def define_limit(pattern_graph):
    num_nodes = pattern_graph.count_nodes()
    num_edges = pattern_graph.count_edges()
    param_matrix = [[0.0 for i in range(num_nodes + 2*num_edges)] for j in range(num_edges)]
    #param_matrix = [[0.0 for i in range(num_nodes)] for j in range(num_edges)]
    b = []                           
    c = [1.0]*num_edges     # DONE
    
    # basic (each ex,  0<=ex<=1)
    #index = node_id + 1
    
    index = 0
    # add for bounds [0,1]
    for edge_id in pattern_graph.get_edges():   # just list of edge_ids
        param_matrix[edge_id][index] = 1.0      # ex <= 1  
        b.append(1.0)
        index += 1
        param_matrix[edge_id][index] = -1.0   # ex >= 0
        b.append(0.0)
        index += 1
    if comment_flag:
        print '==========after adding basis\nmatrix---'
        for each_list in param_matrix:
            print each_list
        print b
        print 'end------'
    
    # fractional
    for node in pattern_graph.get_nodes():   # list of Node instances
        node_id = node.get_node_id()
        node.dump_node()
        # fill it into matrix
        for neighbor in node.get_neighbors():
           edge_id = pattern_graph.get_edge_id(node_id, neighbor)
           print node_id, neighbor, ': ', edge_id
           #param_matrix[edge_id][node_id] = -1.0
           param_matrix[edge_id][index] = -1.0
        b.append(-1.0)
        if comment_flag:   
            print 'matrix---'
            for each_list in param_matrix:
                print each_list
            print b
            print 'end------'
        index+=1

    return param_matrix, b, c
'''
def define_limit_scipy(pattern_graph):
    num_nodes = pattern_graph.count_nodes()
    num_edges = pattern_graph.count_edges()
    param_matrix = [[0.0 for i in range(num_edges)] for j in range(num_nodes)]
    b = []                           
    c = [1.0]*num_edges     # DONE
    
    index = 0
    # fractional
    for node in pattern_graph.get_nodes():   # list of Node instances
        node_id = node.get_node_id()
        #node.dump_node()
        # fill it into matrix
        for neighbor in node.get_neighbors():
           edge_id = pattern_graph.get_edge_id(node_id, neighbor)
           if comment_flag:
               print node_id, neighbor, ': ', edge_id
           param_matrix[index][edge_id] = -1.0
        b.append(-1.0)
        
        if comment_flag:
            print 'matrix---'
            for each_list in param_matrix:
                print each_list
            print b
            print 'end------'
        
        index+=1

    return param_matrix, b, c

def fractional_edge_cover(pattern_graph):
    
    # example for ovopt define matrix
    #A = matrix([ [-1.0, -1.0, 0.0, 1.0], [1.0, -1.0, -1.0, -2.0] ]) # limit(combine with b) -x1+x2, x1+x2, x2, x1-2x2
    #b = matrix([ 1.0, -2.0, 0.0, 4.0 ])  # need to transform to <=, <=1, >=2, >=0, <=4
    #c = matrix([ 2.0, 1.0 ])     # target: 2x1 + x2

    #A,b,c = define_limit(pattern_graph)   # to use ovopt
    A,b,c = define_limit_scipy(pattern_graph)   # to use scipy
    # linear programming
    
    ''' 
    #to use cvopt 
    sol=solvers.lp(matrix(c),matrix(A),matrix(b))
    print(sol['x'])
    count = 0.0
    for value in sol['x']:
        count += value
    print count
    return sol['x']
    '''

    if comment_flag:
        print "start linear programming"
    var_bounds = [(0,1) for i in range(pattern_graph.count_edges())]
    var_bounds = tuple(var_bounds)

    res = linprog(c, A_ub=A, b_ub=b, bounds=var_bounds, options={"disp": True})
    print res

    return res['x']

def split_graph(graph, cover):
    graph_cycles = Graph()
    graph_stars = Graph()

    for index in range(len(cover)):
        src_node_id, target_node_id = graph.translate_edge_id(index) 
        if cover[index] == 0.5:
            graph_cycles.add_edge(src_node_id, target_node_id)
        elif cover[index] == 1.0:
            graph_stars.add_edge(src_node_id, target_node_id)
        elif cover[index] != 0.0:
            print("ALERT: there is edges whose cover values is not in{0, 0.5, 1}")
            exit(1)
        if comment_flag:
            print 'For the cycles',
            graph_cycles.dump_graph()
            print 'For the stars',
            graph_stars.dump_graph()
    
    return graph_cycles, graph_stars

def dfs(graph, start, end):    # graph is instance of Graph, start and end are two node ids
    fringe = [(start, [])]
    while fringe:
        state, path = fringe.pop()
        if path and state == end:
            yield path
            continue
        for next_state in graph.get_neighbors(state):    # iterate neighbors ids of node state
            if next_state in path:
                continue
            fringe.append((next_state, path+[next_state]))

def equal_two_cycles(cycle_1, cycle_2):
    if len(cycle_1) != len(cycle_2):
        return False
    if set(cycle_1) != set(cycle_2):
        return False

    index1 = 0
    index2 = 0
    while cycle_2[index2] != cycle_1[0]:
        index2 += 1

    for index1 in range(len(cycle_1)-1):
        if cycle_1[index1] != cycle_2[index2]:
            return False
        index2 = (index2+1) % (len(cycle_2)-1)   # the last one is meaningless
    
    return True

def find_cycles(graph): 
    # find all cycle paths
    node_ids = []
    for node in graph.get_nodes():
        node_ids.append(node.get_node_id())
    
    if comment_flag:
        print "finding cycles from the graph:"
        graph.dump_graph()

    cycles = [[node]+path  for node in node_ids for path in dfs(graph, node, node)]   # here node is a node id
    if comment_flag:
        print "found", len(cycles), " cycles, they are:"
        print cycles

    # filter cycles: 1) odd 2) different
    uniq_odd_cycles = []
    for cycle in cycles:
        if len(set(cycle))%2 == 1:
            if comment_flag:
                print 'passed odd test', cycle
            flag_uniq = True
            for added_cycle in uniq_odd_cycles:
                if equal_two_cycles(added_cycle, cycle) or equal_two_cycles(added_cycle, list(reversed(cycle))):
                    flag_uniq = False
                    break
            if flag_uniq:
                if comment_flag:
                    print 'passed uniq test'
                uniq_odd_cycles.append(cycle)
        if comment_flag:
            print uniq_odd_cycles
    print("Following are odd, vertex-disjoint cycles within the graph") 
    print uniq_odd_cycles
    # TODO: are these cycles edge-disjoint?
    # if no, choose cycles? cover all edges?
     
    return uniq_odd_cycles

def find_stars(graph):
    cc = graph.connectedComponents() 
    print("Following are stars(connected components)") 
    print(cc) 
 
    return cc

def decompose_graph(pattern_graph, cover):
    # split the pattern to two graphs: 1) all 1s edges and 2) all 1/2 edges
    G_cycles, G_stars = split_graph(pattern_graph, cover)
    
    # find all cycles in the graph with all 1s edges
    cycles = find_cycles(G_cycles)
  
    # find all connected components in the graph with all 1/2 edges
    stars = find_stars(G_stars)

    return cycles, stars    # list of graphs(no cover values)

if __name__=='__main__':
    print '****** This program is to decompose an arbitary pattern graph into odd cycles and stars according to some theory ******'
    if len(sys.argv) < 2:
        print 'to use: python graph_decomp.py input_graph'
        print '        each line of input_graph: source_node target_node'
        exit(1)

    # parse graph
    print "==== parse the input pattern graph ===="
    with open(sys.argv[1], 'r') as graph_in:
        P = parse_graph(graph_in)

    # linear programming to find the Fractional Edge-Cover 
    print "==== finding the fractional edge cover of the graph ===="
    cover = fractional_edge_cover(P)
    print "The optimal fractional cover: ", cover
    for index in range(len(cover)):
        print P.translate_edge_id(index)[0], '--', P.translate_edge_id(index)[1],
        if cover[index] != 0.0:
            if cover[index] == 1.0:
                color = 'red'
            if cover[index] == 0.5:
                color = 'green'
            print '[label=', cover[index], ',color=', color, ']'
        else:
            print ''
    # decompose to odd vertex-disjoint cycles and stars
    print "==== decomposing the support of the fractional edge cover ===="
    cycles, stars = decompose_graph(P, cover)

    with open(sys.argv[1] + '.decomposed', 'w') as graph_out:
        for cycle in cycles:
            graph_out.write(str(cycle) + ' ')
        graph_out.write('\n')
        for star in stars:
            graph_out.write(str(star) + ' ')
