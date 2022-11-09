import sys

with open(sys.argv[1], 'r') as graph_in, open(sys.argv[2], 'w') as graph_out:
	for line in graph_in:
		line_split = line.strip('\n').split(' ')

		graph_out.write(line_split[1] + ' ' + line_split[2] + '\n')
