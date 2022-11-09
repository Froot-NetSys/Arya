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

/* Application 4: counting arbitrary pattern
* sampling times, 0 is enumerating
*/
/*
void pattern_counting(Graph &large_graph, Pattern &graph_pattern, uint64_t sampling_times, bool use_cache) {
    uint64_t lost_edge_cases[100];
    auto start = chrono::high_resolution_clock::now();
    //double estimated_pattern_count = estimate_pattern(large_graph, graph_pattern, sampling_times);
    double this_prob_inverse;
    uint64_t this_sampled_times;
    tie(this_prob_inverse, this_sampled_times) = estimate_pattern(0, large_graph, graph_pattern, sampling_times, lost_edge_cases, use_cache);
    double estimated_pattern_count = this_prob_inverse / this_sampled_times;
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Estimating finds " << estimated_pattern_count << " patterns" << endl;
    cout << "pattern_nums time_consumed(us)\n" << estimated_pattern_count << " " << duration.count() << endl;
    cout << "raw results for merge\n" << this_prob_inverse << " " << this_sampled_times << " " << duration.count()
         << endl;
}
*/

/* Application 5: Multi-threaded counting arbitrary pattern
* sampling times, #threads
* note: sampling times now is #sampling trees, may be larger than sampled profiles
*/
/*
void multi_thread_pattern_counting(Graph &large_graph, Pattern &graph_pattern, uint64_t sampling_times,
                                   uint64_t num_threads, bool use_cache) {
    // overall stats
    vector<Pattern> graph_pattern_vector = {graph_pattern};
    uint64_t lost_edge_cases[100];
    atomic <uint64_t> global_sampled_times(0);                        // global finished sampling times
    // uint64_t sampling_step = min(sampling_times / num_threads, (uint64_t) 10000);   // sampling times for a single step
    uint64_t sampling_step = 10000;
    if (sampling_times < 1000000)
        sampling_step = 1000;
    double global_prob_inverse = 0.0;                          // global prob inverse

    // multi-thread variables
    vector <unique_ptr<thread>> workers(num_threads);      // thread pool
    std::mutex iomutex;                                    // to cleanly print something
    atomic<bool> running_flag(true);                    // signal to stop threads

    for (uint64_t t = 0; t < num_threads; t++) { // Spawn threads
        workers[t] = make_unique<thread>([&, t]() {
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
                tie(this_prob_inverse, this_sampled_times) = estimate_pattern(t, large_graph, graph_pattern_vector,
                                                                              sampling_step, lost_edge_cases, use_cache);

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
*/
/* multi thread sampling, support sampling multiple patterns with the same decomposition basic patterns */
uint64_t total_time = 0;
void multi_thread_multiple_pattern_counting(Graph &large_graph, vector<Pattern> &graph_patterns, vector<uint64_t> &sampling_times_vector,
                                   uint64_t num_threads, bool use_cache) {
    // overall stats
    uint32_t lost_edge_cases[100];
    atomic <uint64_t> global_sampled_times(0);                                                              // global finished sampling times
    vector<atomic <uint64_t>> global_sampled_times_vector(graph_patterns.size());                           // global finished sampling times
    // uint64_t sampling_step = min(sampling_times / num_threads, (uint64_t) 10000);   // sampling times for a single step
    uint64_t sampling_times = *max_element(sampling_times_vector.begin(), sampling_times_vector.end());
    uint64_t sampling_step = 10000;
    if (sampling_times < 1000000)
        sampling_step = 1000;
    vector<double> global_prob_inverse_vector(graph_patterns.size());                          // global prob inverse

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
            vector<double> this_prob_inverse(graph_patterns.size());
            uint64_t this_sampled_times, thread_sampled_times;
            vector<double> thread_prob_inverse(graph_patterns.size());
            thread_sampled_times = 0;

            while (global_sampled_times < sampling_times) {
                tie(this_prob_inverse, this_sampled_times) = estimate_pattern(t, large_graph, graph_patterns,
                                                                              sampling_step, lost_edge_cases, use_cache, this_rand_generator);

                // only change sampled times not prob_inverse, since there is no support of atomic double

                thread_sampled_times += this_sampled_times;
                for (uint32_t k = 0; k < graph_patterns.size(); k++)
                    thread_prob_inverse[k] += global_sampled_times_vector[k] > sampling_times_vector[k] ? 0 : this_prob_inverse[k];

                global_sampled_times += this_sampled_times;
                for (uint32_t k = 0; k < graph_patterns.size(); k++)
                    global_sampled_times_vector[k] += global_sampled_times_vector[k] > sampling_times_vector[k] ? 0 : this_sampled_times;

                // cout << global_sampled_times << " " << sampling_times << endl;
            }

            
                // cout << "=== thread " << thread_id << " is to stop, finished " << thread_sampled_times
                //     << " sampling, prob_inverse: " << thread_prob_inverse << endl;
            for (uint32_t k = 0; k < graph_patterns.size(); k++)
            {
                std::lock_guard <std::mutex> iolock(iomutex);
                global_prob_inverse_vector[k] += thread_prob_inverse[k];
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
    cout << "\n***** main: all threads finished " << endl;
    cout << setprecision(16) << "raw results after merge threads (total_prob_inverse, total_sampled_times): " << endl;
    for (uint32_t k = 0; k < graph_patterns.size(); k++)
        cout << global_prob_inverse_vector[k] << " " << global_sampled_times_vector[k] << endl;

    vector<double> estimated_pattern_count;
    for (uint32_t k = 0; k < graph_patterns.size(); k++)
        estimated_pattern_count.push_back(global_prob_inverse_vector[k] / global_sampled_times_vector[k]);
    cout << "pattern_nums sampling_times\n";
    for (uint32_t k = 0; k < graph_patterns.size(); k++)
        cout << setprecision(16) << estimated_pattern_count[k] << " " << global_sampled_times_vector[k] << endl;
    cout << "total_time_consumed(us) "<< duration.count() << endl;
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
    //uint64_t num_threads = stoi(argv[4]);
    //uint64_t sampling_times = stoi(argv[3]);
    stringstream ss(argv[3]);
    vector <uint64_t> sampling_times_vector;
    while (ss.good()) {
        string substr;
        getline(ss, substr, ',');
        sampling_times_vector.push_back(stoull(substr));
    }
    
    stringstream ss_threads(argv[4]);
    vector <uint64_t> num_threads_vector;
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
    /*
    for (uint32_t k = 0; k < pattern_file_vector.size(); k++) {
	    Pattern graph_pattern(pattern_file_vector[k]);
        graph_pattern.load();
	    // for (uint64_t i = 0; i < sampling_times_vector.size(); i++) {
            for (uint64_t j = 0; j < num_threads_vector.size(); j++) {
                cout << "\n\n**** Current experiment to mine " << pattern_file_vector[k] << " with " << num_threads_vector[j] << " threads, " << sampling_times_vector[k] << " sampling times " << ", use_cache = " << 1 << endl;
                multi_thread_pattern_counting(large_graph, graph_pattern, sampling_times_vector[k], num_threads_vector[j], 1);
        //    }
        }
    }
    */


    // special case: 4-motif
    // 4-motif optimization: calculate 4-motif patterns from one sampling tree
    // 3-star:
    vector<Pattern> graph_patterns;
    Pattern graph_pattern(pattern_file_vector[0]);
    graph_pattern.load();
    graph_patterns.push_back(graph_pattern);
    vector<uint64_t> tmp = {sampling_times_vector[0]};
    // init_cache_index();
    multi_thread_multiple_pattern_counting(large_graph, graph_patterns, tmp, num_threads_vector[0], 0);

/*
    graph_patterns.erase(graph_patterns.begin());
    Pattern graph_pattern1(pattern_file_vector[1]);
    graph_pattern1.load();
    graph_patterns.push_back(graph_pattern1);
    tmp = {sampling_times_vector[1]};
    init_cache_index();
    append_next_cache();
    multi_thread_multiple_pattern_counting(large_graph, graph_patterns, tmp, num_threads_vector[0], 0);
*/
  
    /*
    // 4_chain, 4_cycle, 4_motif_4, 4_motif_5, 4_clique: decomposed as two 1_star
    graph_patterns.erase(graph_patterns.begin());
    sampling_times_vector.erase(sampling_times_vector.begin());
    for (uint32_t k = 1; k < pattern_file_vector.size(); k++) {
        Pattern graph_pattern(pattern_file_vector[k]);
        graph_pattern.load();
        graph_patterns.push_back(graph_pattern);
    }
    multi_thread_multiple_pattern_counting(large_graph, graph_patterns, sampling_times_vector, num_threads_vector[0], 0);
    */
    cout << "\ntotal_running_time = " << total_time << endl;

    return 0;
}
