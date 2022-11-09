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
#include "graph.h"
#include "estimating.h"

using namespace std;


uint32_t total_time = 0;

/* Application 5: Multi-threaded counting arbitrary pattern
* sampling times, #threads
* note: sampling times now is #sampling trees, may be larger than sampled profiles
*/
void multi_thread_pattern_counting(Graph &large_graph, Pattern &graph_pattern, uint64_t sampling_times,
                                   uint32_t num_threads, bool &use_cache_get, bool &use_cache_store) {
    // overall stats
    uint32_t lost_edge_cases[100];
    atomic <uint64_t> global_sampled_times(0);                        // global finished sampling times
    // uint32_t sampling_step = min(sampling_times / num_threads, (uint32_t) 10000);   // sampling times for a single step
    uint32_t sampling_step = 10000;
    if (sampling_times < 1000000)
        sampling_step = 1000;
    double global_prob_inverse = 0.0;                          // global prob inverse

    // multi-thread variables
    vector <unique_ptr<thread>> workers(num_threads);      // thread pool
    std::mutex iomutex;                                    // to cleanly print something
    atomic<bool> running_flag(true);                    // signal to stop threads

    for (uint32_t t = 0; t < num_threads; t++) { // Spawn threads
        workers[t] = make_unique<thread>([&, t]() {
            random_device sd;
            default_random_engine this_rand_generator(sd());
            {
                std::lock_guard <std::mutex> iolock(iomutex);
                // cout << "=== this is thread " << thread_id << " starting to count" << endl;
            }

            // sampling small step by step, until global_sampled_times == sampling_times
            double this_prob_inverse, thread_prob_inverse;
            uint64_t this_sampled_times, thread_sampled_times;
            thread_prob_inverse = 0.0;
            thread_sampled_times = 0;

            while (global_sampled_times < sampling_times) {
                tie(this_prob_inverse, this_sampled_times) = estimate_pattern(t, large_graph, graph_pattern,
                                                                              sampling_step, lost_edge_cases, use_cache_get, use_cache_store, this_rand_generator);

                // only change sampled times not prob_inverse, since there is no support of atomic double
                global_sampled_times += this_sampled_times;

                thread_sampled_times += this_sampled_times;
                thread_prob_inverse += this_prob_inverse;

                // cout << global_sampled_times << " " << sampling_times << endl;
            }

            {
                std::lock_guard <std::mutex> iolock(iomutex);
                // cout << "=== thread " << thread_id << " is to stop, finished " << thread_sampled_times
                //     << " sampling, prob_inverse: " << thread_prob_inverse << endl;
                global_prob_inverse += thread_prob_inverse;
            }
            return;
        });
    }



    //  main thread: test stop and merge results
    auto start = chrono::high_resolution_clock::now();
    for (auto &worker : workers) {
        worker->join();
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    // cout << "***** main: all threads finished " << endl;

    cout << setprecision(16) << "raw results after merge threads (total_prob_inverse, total_sampled_times): " << global_prob_inverse << " "
         << global_sampled_times << endl;
    double estimated_pattern_count = global_prob_inverse / global_sampled_times;
    cout << "pattern_nums time_consumed(us) sampling_times sampling/second\n" << estimated_pattern_count << " " << duration.count()
         << " " << global_sampled_times << " " << global_sampled_times * 1.0 / duration.count() * 1000000 << endl;
    total_time += duration.count();
}

int main(int argc, char *argv[]) {
    // Parse running parameters
    if (argc < 5) {
        cout << "Counting #subgraphs in a large graph \n"
             << "Usage: GraphCounting graph_file subgraph_file(separated by \',\' patterns/4_motif/3_star,patterns/4_motif/4_chain) sampling_times(separated by \',\' 10000,20000) #threads(separated by \',\' 1,4)"
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


    // Create, load graph
    Graph large_graph(graph_file_path);
    large_graph.load();

    // Application: counting arbitrary pattern in the graph

    // multi-threaded:
    bool use_cache_get, use_cache_store;
    for (uint32_t k = 0; k < pattern_file_vector.size(); k++) {
	    Pattern graph_pattern(pattern_file_vector[k]);
        graph_pattern.load();
	    // for (uint32_t i = 0; i < sampling_times_vector.size(); i++) {
            for (uint32_t j = 0; j < num_threads_vector.size(); j++) {
                use_cache_get = use_cache_store = true;
                if (k == 0) 
                    use_cache_get = false;
                if (k == pattern_file_vector.size() -  1)
                    use_cache_store = false;
                cout << "\n\n**** Current experiment to mine " << pattern_file_vector[k] << " with " << num_threads_vector[j] << " threads, " << sampling_times_vector[k] 
                     << " sampling times " << ", use_cache_get = " << use_cache_get << ", use_cache_store = " << use_cache_store << endl;
                if (use_cache_get)
                {
                    init_cache_index();
                    append_next_cache();
                }
                multi_thread_pattern_counting(large_graph, graph_pattern, sampling_times_vector[k], num_threads_vector[j], use_cache_get, use_cache_store);
        //    }
        }
    }
    cout << "total_running_time = " << total_time << endl;

    return 0;
}
