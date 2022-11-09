import os
import sys


check_dir = sys.argv[1]


print 'graph pattern sampling_times estimated_triangles time_consumed(us)'

results = []

for subdir in os.listdir(check_dir):
    #print subdir
    res = ''
    res1 = ''
    res2 = ''
    res3 = ''
    jump = False
    for each_file in os.listdir(check_dir + '/' + subdir):
       if each_file == 'config.json':
           with open(check_dir + '/' + subdir + '/' + each_file, 'r')  as config_file:
               dict_from_file = eval(config_file.read())
               res1 = str(dict_from_file['graph']) + ' ' + str(dict_from_file['pattern']) + ' ' + str(dict_from_file['sampling_times']) + ' ' 
       if each_file == 'running':
           with open(check_dir + '/' + subdir + '/' + each_file, 'r')  as fio_output_file:
               lines = fio_output_file.readlines()
               res2 = lines[-1].strip('\n')

    if jump == True:
        continue
    print res1+res2+res3, subdir
