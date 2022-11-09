#include <iostream>
#include <string.h>
#include <chrono>
#include <atomic>
#include <thread>
#include <pthread.h>
#include <mutex>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iomanip>
#include "mpi.h"
#include <omp.h>
#include "graph.h"
#include "estimating.h"

using namespace std;

#define BUFFER_SIZE 1000

/* Application 6: MPI + OpenMP pattern counting
*/
void MPI_pattern_counting(Graph &large_graph, Pattern &graph_pattern, uint64_t sampling_times,
                                   uint32_t num_threads, uint32_t num_processes) {
    // overall stats
    uint32_t lost_edge_cases[100]; 

    // define sampling_step 
    uint64_t sampling_step = 10000;

    // global variables
    uint64_t global_sampled_times = 0;
    double global_prob_inverse = 0;

    int comm_sz, my_rank;

    MPI_Request recvrqst[BUFFER_SIZE];
    MPI_Status status;
    double recv_prob_inverse[BUFFER_SIZE] = {0};
    uint32_t headptr = 0, tailptr = 0;


    auto start = chrono::high_resolution_clock::now();
#pragma omp parallel num_threads(num_threads)
    {
#pragma omp barrier
#pragma omp master
        {
            int provided;
            MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
            MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
            MPI_Barrier(MPI_COMM_WORLD);
        }
#pragma omp barrier
#pragma omp master  
        {
            while (global_sampled_times < sampling_times)
            {
                if (!my_rank)
                {
                    if (tailptr - headptr < BUFFER_SIZE)
                    {
                        MPI_Irecv(&recv_prob_inverse[tailptr % BUFFER_SIZE], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recvrqst[tailptr % BUFFER_SIZE]);
                        tailptr++;
                    }
                    if (headptr < tailptr)
                    {
                        if (global_sampled_times >= sampling_times) 
                            break;
                        int testflag = 0;
                        MPI_Test(&recvrqst[headptr % BUFFER_SIZE], &testflag, MPI_STATUS_IGNORE); 
                        if (testflag)
                        {
#pragma omp atomic
                            global_prob_inverse += recv_prob_inverse[tailptr % BUFFER_SIZE];
#pragma omp atomic
                            global_sampled_times += sampling_step;
                            headptr++;
                        }
                    }        
                }
                if (my_rank)
                    MPI_Recv(&global_sampled_times, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, &status);
            }
        }
        if (omp_get_thread_num()) {
            MPI_Request sendrqst[BUFFER_SIZE];
            double local_prob_inverse[BUFFER_SIZE] = {0};
            uint64_t local_sampled_times = 0;
            uint32_t headptr = 0, tailptr = 0;

            random_device sd;
            default_random_engine this_rand_generator(sd());

            while (global_sampled_times < sampling_times)
            {
                if (tailptr - headptr < BUFFER_SIZE)
                {
                    tie(local_prob_inverse[tailptr % BUFFER_SIZE], local_sampled_times) = estimate_pattern(large_graph, graph_pattern, sampling_step, lost_edge_cases, this_rand_generator);
                    MPI_Isend(&local_prob_inverse[tailptr % BUFFER_SIZE], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &sendrqst[tailptr % BUFFER_SIZE]);    
                    tailptr++;
                }    
                if (headptr < tailptr)
                {
                    int testflag = 0;
                    MPI_Test(&sendrqst[headptr % BUFFER_SIZE], &testflag, MPI_STATUS_IGNORE);
                    if (testflag) headptr++;
                    if (global_sampled_times >= sampling_times) break;
                }
            }
        }
    }

    
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    if (!my_rank)
    {
        for (uint32_t i = 1; i < num_processes; i++)
            MPI_Send(&global_sampled_times, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
        cout << setprecision(16) << "raw results after merge nodes (total_prob_inverse, total_sampled_times): " << global_prob_inverse << " "
            << global_sampled_times << endl;
        double estimated_pattern_count = global_prob_inverse / global_sampled_times;
        cout << "pattern_nums time_consumed(us) sampling_times sampling/second\n" << estimated_pattern_count << " " << duration.count()
            << " " << global_sampled_times << " " << global_sampled_times * 1.0 / duration.count() * 1000000 << endl;   
    }
    
    MPI_Finalize();
}


int main(int argc, char *argv[]) {
    // Parse running parameters
    if (argc < 6) {
        cout << "Counting #subgraphs in a large graph \n"
             << "Usage: GraphCounting graph_file subgraph_file sampling_times(separated by \',\' 10000,20000) #threads(separated by \',\' 1,4) #processes"
             << endl;
        return 1;
    }

    string graph_file_path = argv[1];
    // string pattern_file_path = argv[2];
    stringstream ss_pattern(argv[2]);
    vector <string> pattern_file_vector;
    while (ss_pattern.good()) {
        string substr;
        getline(ss_pattern, substr, ',');
        pattern_file_vector.push_back(substr);
    }
    //uint32_t num_threads = stoi(argv[4]);
    //uint32_t sampling_times = stoi(argv[3]);
    stringstream ss(argv[3]);
    vector <uint64_t> sampling_times_vector;
    while (ss.good()) {
        string substr;
        getline(ss, substr, ',');
        sampling_times_vector.push_back(stoull(substr));
    }
    
    stringstream ss_threads(argv[4]);
    vector <uint32_t> num_threads_vector;
    while (ss_threads.good()) {
        string substr;
        getline(ss_threads, substr, ',');
        num_threads_vector.push_back(stoi(substr));
    }

    uint32_t num_processes = stoi(argv[5]);

    // Create, load graph
    Graph large_graph(graph_file_path);
    large_graph.load();

    // Application: counting arbitrary pattern in the graph

    // MPI+OpenMP
    for (uint32_t k = 0; k < pattern_file_vector.size(); k++) {
	    Pattern graph_pattern(pattern_file_vector[k]);
        graph_pattern.load();
	    for (uint32_t i = 0; i < sampling_times_vector.size(); i++) {
            for (uint32_t j = 0; j < num_threads_vector.size(); j++) {
                cout << "\n\n**** Current experiment to mine " << pattern_file_vector[k] << " with " << num_threads_vector[j] << " threads, " << sampling_times_vector[i] << " sampling times " << endl;
                MPI_pattern_counting(large_graph, graph_pattern, sampling_times_vector[i], num_threads_vector[j], num_processes);
	        }
        }
    }

    return 0;
}
