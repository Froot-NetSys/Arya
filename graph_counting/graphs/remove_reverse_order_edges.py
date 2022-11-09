import sys

'''
def transform_digraph_to_undrigraph(graph_in, graph_out):
	edges_to_add = {}
	cur_vertex = -1
	cur_vertex_edges = []
	# step 1: load in graph and creat edges that need to be added later
	for line in graph_in:
		line_split = line.strip('\n').split(' ')
		start_vertex = int(line_split[0])
		end_vertex = int(line_split[1])

		# process last finished vertex
		if start_vertex != cur_vertex:
			# join two list, print out cur_vertex_edges
			overall_edges = list(cur_vertex_edges)
			if cur_vertex in edges_to_add:
				overall_edges.extend(x for x in edges_to_add[cur_vertex] if x not in overall_edges)
			
			for another_vertex in overall_edges:
				graph_out.write(str(cur_vertex) + ' ' + str(another_vertex) + '\n')


			cur_vertex = start_vertex		
			cur_vertex_edges = []		


		cur_vertex_edges.append(end_vertex)
		# record edge for later vertex
		if end_vertex not in edges_to_add:
			edges_to_add[end_vertex] = []		
		edges_to_add[end_vertex].append(start_vertex)

	print("after halding every edge")
	for item in sorted(edges_to_add):
		for another_vertex in edges_to_add[item]:
			graph_out.write(str(item) + ' ' + str(another_vertex) + '\n')
'''

def transform_digraph_to_undrigraph(graph_in, graph_out):
	# step 1: load in graph and creat edges that need to be added later
	edges = {}
	for line in graph_in:
		if '#' in line:
			continue
		line_split = line.strip('\n').split('\t')
		# line_split = line.strip('\n').split(' ')
		start_vertex = int(line_split[0])
		end_vertex = int(line_split[1])

		if start_vertex < end_vertex:
			graph_out.write(str(start_vertex) + ' ' + str(end_vertex) + '\n')


with open(sys.argv[1], 'r') as graph_in, open(sys.argv[2], 'w') as graph_out:
	transform_digraph_to_undrigraph(graph_in, graph_out)
