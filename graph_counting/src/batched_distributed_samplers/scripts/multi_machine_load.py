import sys
import subprocess
import threading
    
# read in directory
graph_dir = '' 
num_machines = 0
ips = []
edge_idx = []
graph_path = []

def shcmd(cmd, ignore_error=True):
    print 'Doing:', cmd
    ret = subprocess.call(cmd, shell=True)
    if ignore_error == False and ret != 0:
        print 'Returned', ret, cmd
        exit(ret)
    return ret

class insertThread (threading.Thread):
   def __init__(self, threadID):
      threading.Thread.__init__(self)
      self.threadID = threadID
   
   def run(self):
      i = self.threadID
      print("============== machine " + str(i), graph_path[i], ips[i], edge_idx[i])
        
      # restart memcached
      shcmd("ssh " + ips[i] + " \'sudo service memcached stop ; sudo pkill memcached ; sudo ufw allow 11211 ; memcached -d -m 100000 -p 11211 -u root -M -t 4\' ")

      # sync single_machine_load.py and graph path
      shcmd("scp " + graph_dir + graph_path[i] + " " + ips[i] + ":/users/zz_y/")
      shcmd("scp ./single_machine_load.py " + ips[i] + ":/users/zz_y/")
      # call single_machine_load.py on that node
      # shcmd("ssh " + ips[i] + " \'rm /users/zz_y/" + graph_path[i] + " \' ")
      shcmd("ssh " + ips[i] + " \'cd /users/zz_y/ ; python single_machine_load.py " + graph_path[i] + " " + edge_idx[i][0] + " \' ")
      

if __name__=='__main__':
    # read in directory
    print("python multi_machine_load.py graph_directory(global_subgraph_0..N-1, machine_info)")
    graph_dir = sys.argv[1] + '/'
      
    with open(graph_dir + '/machine_info') as machine_info:
        lines = machine_info.readlines()
        # parse number of machines
        num_machines = int(lines[0].split()[0])
        print("there are " + str(num_machines) + " machines")

        # read in IP and graph path of each machine
        for i in range(0, num_machines):
            line = lines[i+1]
            line_items = line.strip('\n').split()
            
            ips.append(line_items[0])
            edge_idx.append((line_items[3], line_items[4]))
            graph_path.append('global_subgraph_'+str(i))
            # graph_path.append('partition_'+str(i))

            #this_start_v = int(line_items[1])
            #this_end_v = int(line_items[2])

    threads = []
    for i in range(num_machines):
        thread_tmp = insertThread(i)
        threads.append(thread_tmp)
        thread_tmp.start()

    for t in threads:
        t.join()
    print "============ All nodes insert Finished"

