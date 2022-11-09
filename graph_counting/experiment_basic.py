import itertools
from pyreuse.helpers import *
from pyreuse.sysutils.cgroup import Cgroup
import os.path
import time
import json
import os
import datetime
import sys

from pyreuse.sysutils.blocktrace import *
from pyreuse.sysutils.ncq import *
from pyreuse.sysutils.iostat_parser import *

# basis
CL = 64   # size of cache line
KB = 1024
MB = 1024 * KB
GB = 1024 * MB 


# experiment environment



# experiment setup


class Experiment(object):
    def __init__(self):
        # config working directories
        #self.exp_name = '5cycles_citeseer_DFS_test'
        #self.exp_name = '5cycles_100nodes'
        #self.exp_name = 'youtube'
        #self.exp_name = 'citeseer'
        #self.exp_name = 'citeseer_difficult_patterns/7_cycle'
        #self.exp_name = 'citeseer_difficult_patterns/' + sys.argv[1]
        self.exp_name = 'mico_difficult/'
        #self.exp_name = '7cycles_citeseer'
        self.home_dir = '/users/kanwu/single_graph_mining/graph_counting/'
        self.res_dir = self.home_dir + 'results/' + self.exp_name
        self.tmp_dir = '/dev/shm/'
        prepare_dir(self.res_dir)
       
        # tools config
        self.tools_config = {
            'clear_page_cache': True,   # whether clear page cache before each run 
            'blktrace'        : False,   # check block IOs
            'iostat'          : False,  # check ios and cpu/io utilization
            'perf'            : False,  # draw flamegraph
            'sar'             : False   # check page faults
        }

        # experiment config
        config = {
          'graph': ['citeseer.undigraph'],
          'pattern': ['test'],
          'sampling_times': [0, 1000],
        }

        # handle
        self.handle_config(config) 
        print '========================= overall ', len(self.all_configs), ' experiments ============================='
        print '==================== results in ', self.res_dir, ' ============================'
 
    def handle_config(self, config):
        config_dic = list(config.keys())
        config_lists = list(config.values())

        self.all_configs = []
        for element in itertools.product(*config_lists):
            new_config = dict(list(itertools.izip(config_dic, list(element))))
            self.all_configs.append(new_config)
        
    def dump_config(self, config):
        self.cur_exp_dir = self.res_dir + '/' + datetime.now().strftime("%H-%M-%S_%Y%m%d") + sys.argv[1]
        #self.cur_exp_dir = self.res_dir + '/' + datetime.now().strftime("%H-%M-%S_%Y%m%d")
        os.mkdir(self.cur_exp_dir)
        print self.cur_exp_dir
        with open(self.cur_exp_dir + '/config.json', 'w') as config_output:
            json.dump(config, config_output)

    def before_each(self, config):
        print '                ********* Configured with **********'
        print config
        self.dump_config(config)
        # clear page cache
        if self.tools_config['clear_page_cache']:
            shcmd('sync; echo 3 > /proc/sys/vm/drop_caches; sleep 1')
        
 
    def exp(self, config):
        print '              *********** start running ***********'
        cmd = './GraphCounting.out ' + str(config['graph']) + ' ' + str(config['pattern']) + ' ' + str(config['sampling_times']) + ' > ' + self.cur_exp_dir + '/running'

        print cmd
        shcmd(cmd)

    def handle_iostat_out(self, iostat_output):
        print "==== utilization statistics ===="
        stats = parse_batch(iostat_output.read())
        with open(self.cur_exp_dir + '/iostat.out.cpu_parsed', 'w') as parsed_iostat:
            parsed_iostat.write('iowait system user idle \n')
            item_len = average_iowait = average_system = average_user = average_idle = 0
            for item in stats['cpu']:
                parsed_iostat.write(str(item['iowait']) + ' ' + str(item['system']) + ' ' + str(item['user']) + ' ' + str(item['idle']) + '\n')
                if float(item['idle']) > 79:
                    continue
                item_len += 1
                average_iowait += float(item['iowait'])
                average_system += float(item['system'])
                average_user += float(item['user'])
                average_idle += float(item['idle'])
            if item_len > 0:
                print 'iowait  system  user  idle'
                print str(average_iowait/item_len), str(average_system/item_len), str(average_user/item_len), str(average_idle/item_len)
            else:
                print 'seems too idle of CPU'

        with open(self.cur_exp_dir + '/iostat.out.disk_parsed', 'w') as parsed_iostat:
            parsed_iostat.write('r_iops r_bw(MB/s) w_iops w_bw(MB/s) avgrq_sz(KB) avgqu_sz\n')
            item_len = average_rbw = average_wbw = 0
            for item in stats['io']:
                parsed_iostat.write(item['r/s'] + ' ' + item['rMB/s'] + ' ' + item['w/s'] + ' ' + item['wMB/s'] + ' ' + str(float(item['avgrq-sz'])*512/1024) + ' '+ item['avgqu-sz'] +'\n')
                if float(item['rMB/s']) + float(item['wMB/s']) < 20:
                    continue
                item_len += 1
                average_rbw += float(item['rMB/s'])
                average_wbw += float(item['wMB/s'])
            if item_len > 0:
                print str(average_rbw/item_len), str(average_wbw/item_len)
            else:
                print 'seems too idle of Disk'
        print "================================="    

    def after_each(self, config):
        time.sleep(1)
        print '              **************** done ***************'

        # wrapup iostat

    def run(self):
        for config in self.all_configs:
            self.before_each(config)
            self.exp(config)
            self.after_each(config)

if __name__=='__main__':

    configs = []

    # triangles coutning for citeseer
    '''
    config = {
      'graph': ['citeseer.undigraph'],
      'pattern': ['test'],
      'sampling_times': [0, 10, 100, 1000, 2000, 5000, 10000],
    }
    config = {
      'graph': ['mico.undigraph'],
      'pattern': ['test'],
      'sampling_times': [0, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000],
    }
    '''
    
    # 5-cycles counting for citeseer
    '''
    config = {
      'graph': ['citeseer.undigraph'],
      'pattern': ['test_pattern'],
      #'sampling_times': [0],
      'sampling_times': [1000, 10000, 100000, 1000000, 10000000, 20000000, 40000000],
      #'sampling_times': [1000, 10000, 25000, 50000, 75000, 100000, 200000, 500000, 1000000],
    }
    #configs.append(config)

    config = {
      'graph': ['test_100_nodes'],
      'pattern': ['test_pattern'],
      'sampling_times': [1000, 10000, 100000, 1000000, 10000000],
      #'sampling_times': [1000, 10000, 25000, 50000, 75000, 100000, 200000, 500000, 1000000],
    }
    #configs.append(config)
    
    config = {
      'graph': ['./graphs/youtube/youtube.undigraph'],
      'pattern': ['./patterns/5_cycle'],
      #'sampling_times': [1000, 10000, 100000, 1000000, 10000000, 20000000, 40000000],
      'sampling_times': [1000, 10000, 100000, 1000000, 10000000],
    }
    configs.append(config)

    # youtube figure
    config = {
      'graph': ['./graphs/youtube/youtube.undigraph'],
      'pattern': ['./patterns/triangle', './patterns/5_cycle', './patterns/7_cycle', './patterns/triangle_2_star'],
      'sampling_times': [10000, 100000, 1000000, 10000000, 20000000, 40000000, 80000000],
      #'sampling_times': [10000, 100000, 1000000, 10000000],
    }
    configs.append(config)
    '''
    
    # citeseer figure
    config = {
      #'graph': ['./test_graphs/citeseer.undigraph', './test_graphs/mico.undigraph'],
      'graph': ['./test_graphs/mico.undigraph'],
      #'pattern': ['./patterns/triangle', './patterns/5_cycle', './patterns/triangle_2_star', './patterns/triangle_3_star', './patterns/triangle_triangle', './patterns/7_cycle', './patterns/5_cycle_2_star'],
      #'pattern': ['./patterns/2_star', './patterns/3_star', './patterns/3_star_2_star', './patterns/3_star_3_star'],
      #'pattern': ['./patterns/7_cycle','./patterns/5_cycle_2_star', './patterns/triangle_triangle'],
      #'pattern': ['./patterns/3_star_3_star', './patterns/triangle_triangle'],
      'pattern': ['./patterns/triangle_triangle'],
      #'sampling_times': [10000, 100000, 1000000, 10000000, 100000000],
      #'sampling_times': [10000],
      #'sampling_times': [100000000],
      'sampling_times': [50000000],
      'group_id': [sys.argv[1]],
    }
    
    #for i in range(10):
    #    configs.append(config)
    configs.append(config)
    '''

    # mico figure
    config = {
      'graph': ['./test_graphs/mico.undigraph'],
      'pattern': ['./patterns/triangle', './patterns/5_cycle', './patterns/triangle_triangle'],
      'sampling_times': [10000, 100000, 1000000, 10000000, 20000000, 40000000, 80000000],
      #'sampling_times': [160000000, 320000000, 640000000],
    }
    configs.append(config)
    configs.append(config)
    configs.append(config)
    configs.append(config)
    configs.append(config)
    '''


    # run each config
    for each_config in configs:
        exp = Experiment()
        exp.handle_config(each_config)
        exp.run()
