import sys
import threading
import pylibmc

n_threads = 4

edges = []
edges_len = 0
start_edge_id = 0

class insertThread (threading.Thread):
   def __init__(self, threadID):
      threading.Thread.__init__(self)
      self.threadID = threadID
   
   def run(self):
      print "Starting thread " + str(self.threadID)
      # creat a client to connect to local memcached
      mc = pylibmc.Client(["127.0.0.1"], binary=True, behaviors={"tcp_nodelay": True, "ketama": True})
      
      idx = self.threadID
      local_count = 0
      while idx < edges_len:
          vertexes = edges[idx].strip('\n').split()
          mc.set("edge-"+str(start_edge_id + idx), vertexes[0] + " " + vertexes[1], 1000000)


          # test insertion
          # print("inserted " + str(start_edge_id + idx))
          '''
          value = mc.get("edge-" + str(idx))
          # print(value.split()[0], value.split()[1])
          assert value.split()[0] == vertexes[0] and value.split()[1] == vertexes[1]
          '''

          idx += n_threads
          local_count += 1
          if local_count % 100000 == 0:
              print("thread " + str(self.threadID) + " finished " + str(local_count))

      return

if __name__=='__main__':
    print("python single_machine_load.py graph_file_path start_edge_id(0 if this is the only node)")
    graph_path = sys.argv[1]
    start_edge_id = int(sys.argv[2])

    print("insert " + graph_path + " edge id starting from " + str(start_edge_id))
    
    mc = pylibmc.Client(["127.0.0.1"], binary=True, behaviors={"tcp_nodelay": True, "ketama": True})
    mc.flush_all()

    # read in edges
    with open(graph_path) as graph:
        edges = graph.readlines()
        edges_len = len(edges)
  
    # Step 1: multi-threading insert edges
    threads = []
    for i in range(n_threads):
        thread_tmp = insertThread(i)
        threads.append(thread_tmp)
        thread_tmp.start()

    for t in threads:
        t.join()
    print "Insert Edges Finished"

    # Step 2: insert vertex degrees
    cur_v = sys.maxint
    cur_v_neighbors = []
    for edge in edges:
        v_start = edge.strip('\n').split()[0]
        v_end = edge.strip('\n').split()[1]

        if v_start != cur_v:
            if cur_v != sys.maxint:
                # insert adj list to memcached
                adj_list_str = ' '.join(cur_v_neighbors)
                mc.set("v-"+str(cur_v), adj_list_str, 1000000)

                # test edges
                # val = mc.get('v-'+str(cur_v))
                # print('v-'+str(cur_v) + " : " + val)

            cur_v = v_start
            cur_v_neighbors = []

        
        cur_v_neighbors.append(v_end)


    adj_list_str = ' '.join(cur_v_neighbors)
    mc.set("v-"+str(cur_v), adj_list_str, 1000000)
    # test edges
    # val = mc.get('v-'+str(cur_v))
    # print('v-'+str(cur_v) + " : " + val)
    
    print "Insert Vertexes Adj List Finished"
