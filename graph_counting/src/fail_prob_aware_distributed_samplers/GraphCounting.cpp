#include <iostream>
#include <iomanip>
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
#include <vector>
#include <string>
#include "graph.h"
#include "estimating.h"


using namespace std;

/* Application 5: Multi-threaded counting arbitrary pattern
* sampling times, #threads
* note: sampling times now is #sampling trees, may be larger than sampled profiles
*/
double multi_thread_pattern_counting(Graph &large_graph, Pattern &graph_pattern, uint64_t sampling_times,
                                     uint32_t num_threads, string machine_info_file_path) {
    // overall stats
    atomic<uint64_t> global_cross_sampled_times(0);                        // global finished sampling times
    atomic<uint64_t> global_intra_sampled_times(0);
    double global_cross_prob_inverse = 0.0;                          // global prob inverse
    double global_intra_prob_inverse = 0.0;

    uint32_t sampling_step = sampling_times / num_threads;
    vector<uint32_t> worker_sampling_times(num_threads, sampling_step);
    for (uint32_t t = 0; t < sampling_times % num_threads; t++)
        worker_sampling_times[t] += 1;

    // multi-thread variables
    vector<unique_ptr<thread>> workers(num_threads);      // thread pool
    std::mutex iomutex;                                    // to cleanly print something

    for (uint32_t t = 0; t < num_threads; t++) { // Spawn threads
        workers[t] = make_unique<thread>([&, t]() {
            uint32_t thread_id = t;
            random_device sd;
            default_random_engine this_rand_generator(sd());
            {
                std::lock_guard<std::mutex> iolock(iomutex);
                cout << "=== this is thread " << thread_id << " starting to count" << endl;
                // create libmemcached client for each thread
                // use sockets per memcached
                large_graph.load_machine_info(machine_info_file_path, thread_id);
            }

            // sampling small step by step, until global_sampled_times == sampling_times
            double cross_prob_inverse, intra_prob_inverse;
            uint64_t num_cross_partition, num_intra_partition;

            tie(cross_prob_inverse, num_cross_partition, intra_prob_inverse, num_intra_partition) = estimate_pattern(
                    thread_id, large_graph, graph_pattern, worker_sampling_times[t], this_rand_generator);
            // tie(cross_prob_inverse, num_cross_partition, intra_prob_inverse, num_intra_partition) = estimate_stars(thread_id, large_graph, graph_pattern, worker_sampling_times[t]);

            // only change sampled times not prob_inverse, since there is no support of atomic double
            // global_sampled_times += worker_sampling_times[t];
            global_cross_sampled_times += num_cross_partition;
            global_intra_sampled_times += num_intra_partition;

            // cout << global_sampled_times << " " << sampling_times << endl;

            {
                std::lock_guard<std::mutex> iolock(iomutex);
                /*
                 cout << "=== thread " << thread_id << " is to stop, finished " << thread_sampled_times
                     << " sampling, prob_inverse: " << thread_prob_inverse << endl;
                */
                global_cross_prob_inverse += cross_prob_inverse;
                global_intra_prob_inverse += intra_prob_inverse;
            }

            return;
        });
    }



    //  main thread: test stop and merge results
    auto start = chrono::high_resolution_clock::now();
    for (auto &worker: workers) {
        worker->join();
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);

    // cout << "global_cross_prob_inverse global_cross_sampled_times global_intra_prob_inverse global_intra_sampled_times" << endl;
    //cout << global_cross_prob_inverse << " " << global_cross_sampled_times << " " << global_intra_prob_inverse << " " << global_intra_sampled_times << " " << endl;

    cout << "cross_partition_count intra_partition_count" << endl;
    cout << setprecision(16) << global_cross_prob_inverse / sampling_times << " "
         << global_intra_prob_inverse / sampling_times << endl;
    cout << "total_pattern_count estimated_count_ratio time_consumed" << endl;
    cout << (global_cross_prob_inverse + global_intra_prob_inverse) / sampling_times << " "
         << (global_cross_prob_inverse + global_intra_prob_inverse) / global_intra_prob_inverse << " "
         << duration.count() << endl;

    return (global_cross_prob_inverse + global_intra_prob_inverse) / global_intra_prob_inverse;



    // cout << "***** main: all threads finished " << endl;

    // cout << "raw results after merge threads (total_prob_inverse, total_sampled_times): " << global_prob_inverse << " "
    //      << global_sampled_times << endl;
    // double estimated_pattern_count = global_prob_inverse / global_sampled_times;
    // cout << "pattern_nums time_consumed(us) sampling_times sampling/second\n" << estimated_pattern_count << " " << duration.count()
    //      << " " << global_sampled_times << " " << global_sampled_times * 1.0 / duration.count() * 1000000 << endl;
    // return make_tuple(estimated_pattern_count, duration.count());
}


double multi_thread_batch_counting(Pattern &graph_pattern, uint64_t sampling_times, uint32_t num_threads,
                                   string machine_info_file_path) {
    // overall stats
    atomic<uint64_t> global_cross_sampled_times(0);                        // global finished sampling times
    atomic<uint64_t> global_intra_sampled_times(0);
    double global_cross_prob_inverse = 0.0;                          // global prob inverse
    double global_intra_prob_inverse = 0.0;

    uint32_t sampling_batch = 50;
    uint32_t per_thread_sampling_times = sampling_times / num_threads + 1;

    // multi-thread variables
    vector<unique_ptr<thread>> workers(num_threads);      // thread pool
    std::mutex iomutex;                                    // to cleanly print something

    for (uint32_t t = 0; t < num_threads; t++) { // Spawn threads
        workers[t] = make_unique<thread>([&, t]() {
            uint32_t thread_id = t;
            Graph large_graph;

            uint32_t local_sampled_times = 0;
            uint64_t local_cross_sampled_times = 0;
            uint64_t local_intra_sampled_times = 0;
            double local_cross_prob_inverse = 0;
            double local_intra_prob_inverse = 0;

            double cross_prob_inverse, intra_prob_inverse;
            uint64_t num_cross_partition, num_intra_partition;

            random_device sd;
            default_random_engine this_rand_generator(sd());
            
            // Step 1: create connections to each memcached node
            {
                std::lock_guard<std::mutex> iolock(iomutex);
                cout << "=== this is thread " << thread_id << " starting to count" << endl;
                // create sockets for each memcached node
                large_graph.load_machine_info(machine_info_file_path, 0);
            }

            // sampling small step by step, until global_sampled_times == sampling_times
            while (local_sampled_times < per_thread_sampling_times) {
                // tie(cross_prob_inverse, num_cross_partition, intra_prob_inverse, num_intra_partition) = estimate_pattern(thread_id, large_graph, graph_pattern, worker_sampling_times[t]);
                // tie(cross_prob_inverse, num_cross_partition, intra_prob_inverse, num_intra_partition) = estimate_stars(0,large_graph,graph_pattern, sampling_batch);
                // tie(cross_prob_inverse, num_cross_partition, intra_prob_inverse, num_intra_partition) = estimate_odd_cycles(0,large_graph,graph_pattern, sampling_batch);
                tie(cross_prob_inverse, num_cross_partition, intra_prob_inverse,
                    num_intra_partition) = estimate_pattern_batch(0, large_graph, graph_pattern, sampling_batch, this_rand_generator);

                local_cross_sampled_times += num_cross_partition;
                local_intra_sampled_times += num_intra_partition;
                local_cross_prob_inverse += cross_prob_inverse;
                local_intra_prob_inverse += intra_prob_inverse;

                local_sampled_times += sampling_batch;
            }

            {
                std::lock_guard<std::mutex> iolock(iomutex);
                global_cross_sampled_times += local_cross_sampled_times;
                global_intra_sampled_times += local_intra_sampled_times;
                global_cross_prob_inverse += local_cross_prob_inverse;
                global_intra_prob_inverse += local_intra_prob_inverse;
            }


            return;
        });
    }



    //  main thread: test stop and merge results
    auto start = chrono::high_resolution_clock::now();
    for (auto &worker: workers) {
        worker->join();
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);

    cout << "cross_partition_count intra_partition_count" << endl;
    cout << setprecision(16) << global_cross_prob_inverse / sampling_times << " "
         << global_intra_prob_inverse / sampling_times << endl;
    cout << "total_pattern_count estimated_count_ratio time_consumed" << endl;
    cout << (global_cross_prob_inverse + global_intra_prob_inverse) / sampling_times << " "
         << (global_cross_prob_inverse + global_intra_prob_inverse) / global_intra_prob_inverse << " "
         << duration.count() << endl;

    return (global_cross_prob_inverse + global_intra_prob_inverse) / global_intra_prob_inverse;
}

double batch_counting(Pattern &graph_pattern, uint64_t sampling_times, uint32_t num_threads,
                      string machine_info_file_path) {
    // overall stats
    auto start = chrono::high_resolution_clock::now();

    atomic<uint64_t> global_cross_sampled_times(0);                        // global finished sampling times
    atomic<uint64_t> num_succeed(0);
    double global_cross_prob_inverse = 0.0;                          // global prob inverse
    double global_intra_prob_inverse = 0.0;

    uint32_t local_sampled_times = 1;

    random_device sd;
    default_random_engine this_rand_generator(sd());

    // Step 1: create connections to each memcached node
    Graph large_graph;
    large_graph.load_machine_info(machine_info_file_path, 0);

    // sampling small step by step, until global_sampled_times == sampling_times
    double cross_prob_inverse, intra_prob_inverse;
    uint64_t num_cross_partition, num_intra_partition;
    bool flag = 0;

    while (local_sampled_times < sampling_times) {
        tie(cross_prob_inverse, num_cross_partition, intra_prob_inverse, num_intra_partition) = estimate_pattern_batch(
                0, large_graph, graph_pattern, sampling_batch, this_rand_generator);

        global_cross_sampled_times += num_cross_partition;
        global_cross_prob_inverse += cross_prob_inverse;
        num_succeed += intra_prob_inverse;

        local_sampled_times += sampling_batch;
        if (!flag && local_sampled_times > sampling_times * 10)
        {
            adjust_sampling_block_order(graph_pattern);
            flag = 1;
        }
    }

    //  main thread: test stop and merge results
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);

    cout << "total_pattern_count total_counter sampling_times time_consumed succeed_times fail_probability" << endl;
    cout << setprecision(16)
         << (global_cross_prob_inverse + global_intra_prob_inverse) / global_cross_sampled_times
         << " " << (global_cross_prob_inverse + global_intra_prob_inverse)
         << " " << global_cross_sampled_times
         << " " << duration.count()
         << " " << num_succeed
         << " " << 1.0 * num_succeed / global_cross_sampled_times
         << endl;



    return 0.0;
}

int main(int argc, char *argv[]) {
    // Parse running parameters
    if (argc < 6) {
        cout << "Counting #subgraphs in a large graph \n"
             << "Usage: GraphCounting graph_file(not used) subgraph_file machine_info_file sampling_times #threads "
             << endl;
        return 1;
    }

    string graph_file_path = argv[1];
    string pattern_file_path = argv[2];
    string machine_info_file_path = argv[3];
    uint32_t sampling_times = stoi(argv[4]);
    uint32_t num_threads = stoi(argv[5]);


    // Application: counting arbitrary pattern in the graph
    Pattern graph_pattern(pattern_file_path);
    graph_pattern.load();

    // multi-threaded memcached global sampling
    double estimated_count_ratio = 0; // total / intra #pattern
    // Graph large_graph;
    // estimated_count_ratio = multi_thread_pattern_counting(large_graph, graph_pattern, sampling_times, num_threads, machine_info_file_path);
    // estimated_count_ratio = multi_thread_batch_counting(graph_pattern, sampling_times, num_threads,
    //                                                    machine_info_file_path);
    estimated_count_ratio = batch_counting(graph_pattern, sampling_times, num_threads,
                                                        machine_info_file_path);

    return 0;
}
