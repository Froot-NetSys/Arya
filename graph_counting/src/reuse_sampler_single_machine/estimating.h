// estimating.h
#ifndef ESTIMATING_H
#define ESTIMATING_H

#include <iostream>
#include "graph.h"
#include <cmath>
#include <set>
#include <queue>
#include <map>
#include <unordered_map>
#include <mutex>
#include "assert.h"

#define MAX_THREAD_NUM 40

using namespace std;

/* SamplingTree strcuture
 * tree of nodes (vertex, probability)
*/
struct SamplingTreeNode {
    SamplingTreeNode(NodeType type_, vector <uint32_t> &vertexes_, double probability_) : node_type(type_),
                                                                                          probability(probability_) {
        for (uint32_t i = 0; i < vertexes_.size(); i++) {
            vertexes.push_back(vertexes_[i]);
        }
    }

    void print() {
        cout << "Node(" << NodeTypeStr[node_type] << ", " << children.size() << ", " << probability << ", [ ";
        for (auto it = vertexes.begin(); it != vertexes.end(); it++)
            cout << *it << " ";
        cout << "]) ";
    }

    NodeType node_type;   // CYCLE_ROOT, CYCLE_LEAF, STAR_ROOT, STAR_LEAF
    vector <uint32_t> vertexes;
    vector<struct SamplingTreeNode *> children;
    double probability;
};


/* Cache (reuse) of k-cycle and k-star sampling tree nodes
   per-thread cache to avoid lock
*/
unordered_map<uint32_t, vector<struct SamplingTreeNode *>> cycle_cache[MAX_THREAD_NUM];
unordered_map<uint32_t, vector<struct SamplingTreeNode *>> star_cache[MAX_THREAD_NUM];
unordered_map<uint32_t, uint64_t> cycle_cache_index[MAX_THREAD_NUM];
unordered_map<uint32_t, uint64_t> star_cache_index[MAX_THREAD_NUM];
unordered_map<uint32_t, vector<struct SamplingTreeNode *>> cycle_cache_next[MAX_THREAD_NUM];
unordered_map<uint32_t, vector<struct SamplingTreeNode *>> star_cache_next[MAX_THREAD_NUM];

/* copy node data from root_b to root_a */
struct SamplingTreeNode * copy_sampling_tree_node(struct SamplingTreeNode * root_b) {
    struct SamplingTreeNode * root_a = new SamplingTreeNode(root_b->node_type,
                                                            root_b->vertexes,
                                                            root_b->probability);
    assert(root_a != nullptr);
    for (uint32_t j = 0; j < root_b->children.size(); j++) {
        root_a->children.push_back(new SamplingTreeNode(root_b->children[j]->node_type, 
                                                        root_b->children[j]->vertexes, 
                                                        root_b->children[j]->probability));
    }
    return root_a;
}

void store_to_cache(unordered_map<uint32_t, vector<struct SamplingTreeNode *>> * cycle_cache, 
    unordered_map<uint32_t, vector<struct SamplingTreeNode *>> * star_cache,
    uint32_t thread_id, struct SamplingTreeNode * root, uint32_t basic_pattern_para) {
        struct SamplingTreeNode * root_cache = copy_sampling_tree_node(root);
        if (root->node_type == StarRoot) {
            // std::lock_guard <std::mutex> lock(write_mutex);
            vector<struct SamplingTreeNode *> &vec = star_cache[thread_id][basic_pattern_para];
            vec.push_back(root_cache);
            // cout << star_cache[basic_pattern_para].size() << endl;
        }
        else if (root->node_type == CycleRoot) {
            //std::lock_guard <std::mutex> lock(write_mutex);
            vector<struct SamplingTreeNode *> &vec = cycle_cache[thread_id][basic_pattern_para];
            vec.push_back(root_cache);
            // cout << cycle_cache[basic_pattern_para].size() << endl;
        }
}

struct SamplingTreeNode * get_from_cache(unordered_map<uint32_t, vector<struct SamplingTreeNode *>> * cycle_cache, 
    unordered_map<uint32_t, vector<struct SamplingTreeNode *>> * star_cache,
    uint32_t thread_id, NodeType type, uint32_t basic_pattern_para) {
        if (type == StarRoot && star_cache[thread_id].find(basic_pattern_para) != star_cache[thread_id].end()) {
            if (star_cache[thread_id].at(basic_pattern_para).size() == 0)
                return nullptr;
            if (star_cache_index[thread_id].find(basic_pattern_para) == star_cache_index[thread_id].end())
                star_cache_index[thread_id][basic_pattern_para] = 0;
            uint64_t & index = star_cache_index[thread_id][basic_pattern_para];
            if (index >= star_cache[thread_id].at(basic_pattern_para).size())
                return nullptr;
            return star_cache[thread_id].at(basic_pattern_para)[index++]; 
        }
        else if (type == CycleRoot && cycle_cache[thread_id].find(basic_pattern_para) != cycle_cache[thread_id].end()) {
            if (cycle_cache[thread_id].at(basic_pattern_para).size() == 0)
                return nullptr;
            if (cycle_cache_index[thread_id].find(basic_pattern_para) == cycle_cache_index[thread_id].end())
                cycle_cache_index[thread_id][basic_pattern_para] = 0;
            uint64_t & index = cycle_cache_index[thread_id][basic_pattern_para];
            if (index >= cycle_cache[thread_id].at(basic_pattern_para).size())
                return nullptr;
            return cycle_cache[thread_id].at(basic_pattern_para)[index++];
        }
    return nullptr;
}


/* Sample a odd-cycle
* from graph, <u1, u2, ..., w1>
* node([(v1, v2), (v3, v4), ...]) -> k x w1
*/
SamplingTreeNode *odd_cycle_sampler(uint32_t thread_id, Graph &graph, uint32_t cycle_length, bool use_cache, default_random_engine & rand_generator) {
    SamplingTreeNode *root = nullptr;
    if (use_cache) {
        SamplingTreeNode * root_cache = nullptr;
        root_cache = get_from_cache(cycle_cache, star_cache, thread_id, CycleRoot, cycle_length);
        // TODO: copy to new root
        if (root_cache != nullptr) {
            root = copy_sampling_tree_node(root_cache);
        }
    }
    if (!use_cache || root == nullptr) {
        //cout << "====== Start to sample a " << cycle_length << "-cycles" << endl;
        // default_random_engine rand_generator(std::chrono::system_clock::now().time_since_epoch().count());
        uint32_t ue_id, d_ue;

        // step 1: sample cycle_length/2 edges (u1, v1), (x, x), (x, x), ... u1 < v1
        vector <uint32_t> sample_edge_vector;
        Edge e1 = graph.edge_sampling(rand_generator);
        if (graph.vertex_before(e1.v_start, e1.v_end)) {
            sample_edge_vector.push_back(e1.v_start);
            sample_edge_vector.push_back(e1.v_end);
        } else {
            sample_edge_vector.push_back(e1.v_end);
            sample_edge_vector.push_back(e1.v_start);
        }

        for (uint32_t j = 0 + 1; j < cycle_length / 2; j++) {
            Edge sampled_edge = graph.edge_sampling(rand_generator);
            sample_edge_vector.push_back(sampled_edge.v_start);
            sample_edge_vector.push_back(sampled_edge.v_end);
        }

        root = new SamplingTreeNode(CycleRoot,
                                    sample_edge_vector,
                                    pow(graph.total_edges_, (cycle_length / 2)) / 2.0);
        assert(root != nullptr);


        // step 2: sample w1, a neighbor of u1
        ue_id = sample_edge_vector[0];
        d_ue = graph.degree(ue_id);

        // cout << "cycle u1 degree: " << d_ue << endl;

        uniform_int_distribution <uint32_t> next_edge_distribution(0, d_ue - 1);

        for (uint32_t j = 0; j < ceil(d_ue / graph.sqrt_m); j++) {
            // step 2: sample an vertex from u1's neighbor
            vector <uint32_t> sampled_vertex_id;
            sampled_vertex_id.push_back(graph.retrieve_neighbor(ue_id, next_edge_distribution(rand_generator)));
            root->children.push_back(new SamplingTreeNode(CycleLeave, sampled_vertex_id, d_ue));
        }
        if (use_cache)
            store_to_cache(cycle_cache_next, star_cache_next, thread_id, root, cycle_length);
    }
    return root;
}

/* Sample a star
* from graph, <center, w1, w2, ...>
*/
SamplingTreeNode *star_sampler(uint32_t thread_id, Graph &graph, uint32_t petal_num, bool use_cache, default_random_engine & rand_generator) {
    SamplingTreeNode *root = nullptr;
    if (use_cache) {
        SamplingTreeNode * root_cache = nullptr;
        root_cache = get_from_cache(cycle_cache, star_cache, thread_id, StarRoot, petal_num);
        if (root_cache != nullptr) {
            root = copy_sampling_tree_node(root_cache);
        }
    }
    if (!use_cache || root == nullptr) {
        //cout << "====== Start to sample an " << petal_num << "-star" << endl;
        // default_random_engine rand_generator(std::chrono::system_clock::now().time_since_epoch().count());

        // Step 1: choose a vertex, proportional to its degree, choose an edge, its first vertex
        Edge e1 = graph.edge_sampling(rand_generator);
        uint32_t center_vertex_id = e1.v_start;
        uint32_t center_vertex_degree = graph.degree(center_vertex_id);
        vector <uint32_t> center_vertex;
        center_vertex.push_back(center_vertex_id);
        root = new SamplingTreeNode(StarRoot,
                                    center_vertex,
                                    1.0 * graph.total_edges_ / center_vertex_degree);

        // Step 2: sample petal_num vertexes from N(X), without replacement
        set <uint32_t> vertex_set;
        uniform_int_distribution <uint32_t> next_edge_distribution(0, center_vertex_degree - 1);

        vector <uint32_t> sampled_petals;
        if (center_vertex_degree >= petal_num) {
            // indicate can find enough petals
            for (uint32_t i = 0; i < petal_num; i++) {
                uint32_t rand_index = next_edge_distribution(rand_generator);
                while (vertex_set.find(rand_index) != vertex_set.end()) {
                    rand_index = next_edge_distribution(rand_generator);
                }
                vertex_set.insert(rand_index);
                uint32_t sampled_vertex_id = graph.retrieve_neighbor(center_vertex_id, rand_index);
                sampled_petals.push_back(sampled_vertex_id);
            }
        }

        // cout << "star vertex degree: " << center_vertex_degree << " , petal num: " << petal_num << endl;
        uint64_t this_count = 1;
        if (petal_num < center_vertex_degree) {
            for (uint32_t i = center_vertex_degree; i > center_vertex_degree - petal_num; i--)
                this_count = this_count * i;

            for (uint32_t i = 2; i <= petal_num; i++)
                this_count = this_count / i;
        }

        root->children.push_back(new SamplingTreeNode(StarLeave, sampled_petals, 1.0 * this_count));

        // cout << center_vertex_id << " " << center_vertex_degree
        //        << " " << root->probability << " " <<  root->children[0]->probability
        //        << " " << root->probability * root->children[0]->probability << endl;
        if (use_cache)
            store_to_cache(cycle_cache_next, star_cache_next, thread_id, root, petal_num);
    }
    return root;
}

/* BFS Iterate sampling trees to calculate a value
*/
void printTree(SamplingTreeNode *root) {
    // data structure
    queue < SamplingTreeNode * > to_visit_list;
    int depth = 0;

    // insert root
    to_visit_list.push(root);

    // BFS to print
    while (to_visit_list.size() > 0) {
        cout << "level " << depth << ": ";
        int sz = to_visit_list.size();
        for (int i = 0; i < sz; i++) {
            // pop the front
            SamplingTreeNode *next_print = to_visit_list.front();
            if (next_print == nullptr) {
                cout << "!!!!!! get a nullptr" << endl;
            }
            next_print->print();

            // expand next level
            for (uint32_t j = 0; j < next_print->children.size(); j++) {
                to_visit_list.push(next_print->children[j]);
            }

            to_visit_list.pop();
        }
        cout << endl;
        depth++;
    }
    return;
}

/* Create sampling trees
*/
SamplingTreeNode *construct_sampling_tree(uint32_t thread_id, Graph &graph, Pattern &pattern, bool use_cache, default_random_engine & rand_generator) {
    SamplingTreeNode *root = nullptr;
    vector < SamplingTreeNode * > leaves;
    vector < SamplingTreeNode * > new_leaves;

    // Step 1: iterate odd-cycles to add 2-levels
    // cout << "grow tree with cycles" << endl;
    for (uint32_t i = 0; i < pattern.cycle_num; i++) {
        if (i == 0) {
            root = odd_cycle_sampler(thread_id, graph, pattern.cycles[i].size(), use_cache, rand_generator);
            for (uint32_t j = 0; j < root->children.size(); j++) {
                leaves.push_back(root->children[j]);
            }
        } else {
            // extend every leaf by a cycle tree
            for (uint32_t j = 0; j < leaves.size(); j++) {
                leaves[j]->children.push_back(odd_cycle_sampler(thread_id, graph, pattern.cycles[i].size(), use_cache, rand_generator));
            }

            // update leaves and new_leaves
            new_leaves.clear();
            for (uint32_t j = 0; j < leaves.size(); j++) {
                SamplingTreeNode *last_level_root = leaves[j]->children[0];
                for (uint32_t k = 0; k < last_level_root->children.size(); k++) {
                    new_leaves.push_back(last_level_root->children[k]);
                }
            }
            leaves = new_leaves;
        }
    }

    // Step 2: iterate stars to add 2-levels
    // cout << "grow tree with stars" << endl;
    for (uint32_t i = 0; i < pattern.star_num; i++) {
        if (root == nullptr) {
            root = star_sampler(thread_id, graph, pattern.stars[i].size() - 1, use_cache, rand_generator);
            for (uint32_t j = 0; j < root->children.size(); j++) {
                leaves.push_back(root->children[j]);
            }
        } else {
            // extend every leaf by a cycle tree
            for (uint32_t j = 0; j < leaves.size(); j++) {
                leaves[j]->children.push_back(star_sampler(thread_id, graph, pattern.stars[i].size() - 1, use_cache, rand_generator));
            }

            // update leaves and new_leaves
            new_leaves.clear();
            for (uint32_t j = 0; j < leaves.size(); j++) {
                SamplingTreeNode *last_level_root = leaves[j]->children[0];
                for (uint32_t k = 0; k < last_level_root->children.size(); k++) {
                    new_leaves.push_back(last_level_root->children[k]);
                }
            }
            leaves = new_leaves;
        }
    }

    return root;
}

/* Testing whether profiles -> actual copy of patterns
* count probabilities as well
*/
bool test_circle_copy(Graph &graph, vector <uint32_t> &e_vertexes, uint32_t w) {
    // test u1 is the smallest vertex
    uint32_t u1 = e_vertexes[0];
    uint32_t v1 = e_vertexes[1];
    for (uint32_t i = 1; i < e_vertexes.size(); i++) {
        if (!graph.vertex_before(u1, e_vertexes[i]))
            return false;
    }

    if (!graph.vertex_before(u1, w))
        return false;

    // test v1 < w
    if (!graph.vertex_before(v1, w))
        return false;

    // test edges of (v1,u2), ... (vk, w)
    if (!graph.neighbor_test(e_vertexes.back(), w))
        return false;

    for (uint32_t i = 1; i < e_vertexes.size() - 1; i += 2) {
        if (!graph.neighbor_test(e_vertexes[i], e_vertexes[i + 1]))
            return false;
    }

    return true;
}

/* Major helper function to test a sampled path
*/
void DFS_judge_remaining_edges(uint32_t idx, Graph &graph, Pattern &pattern,
                               vector <vector<uint32_t>> &cycle_vertexes, vector <vector<uint32_t>> &star_vertexes,
                               map <uint32_t, uint32_t> &v_map, uint64_t &found_instances, uint32_t *lost_edge_cases) {
    // corner cases
    if (idx == cycle_vertexes.size() + star_vertexes.size()) {

        // judge remaining edges
        bool test_succeed = true;
        uint32_t lost_edges = 0;
        for (uint32_t i = 0; i < pattern.other_edges_.size(); i++) {
            uint32_t test_v_start = v_map.at(pattern.other_edges_[i].v_start);
            uint32_t test_v_end = v_map.at(pattern.other_edges_[i].v_end);

            if (!graph.neighbor_test(test_v_start, test_v_end)) {
                test_succeed = false;
                lost_edges += 1;
                break; // early return 
            }
        }

        if (test_succeed) {
            found_instances += 1;
        }
        lost_edge_cases[lost_edges] += 1;
        return;
    }

    // normal grow mapping
    if (idx < cycle_vertexes.size()) {
        for (uint32_t i = 0; i < pattern.candidate_cycle_v_maps[idx].size(); i++) {
            // map candidate_cycle_v_maps[idx][i] to cycle_vertexes[idx]
            for (uint32_t j = 0; j < pattern.candidate_cycle_v_maps[idx][i].size(); j++) {
                v_map.emplace(pattern.candidate_cycle_v_maps[idx][i][j], cycle_vertexes[idx][j]);
            }

            // DFS_judge_further
            DFS_judge_remaining_edges(idx + 1, graph, pattern, cycle_vertexes, star_vertexes, v_map, found_instances, lost_edge_cases);

            // restore
            for (uint32_t j = 0; j < pattern.candidate_cycle_v_maps[idx][i].size(); j++) {
                v_map.erase(pattern.candidate_cycle_v_maps[idx][i][j]);
            }
        }
    } else {
        uint32_t star_idx = idx - cycle_vertexes.size();
        for (uint32_t i = 0; i < pattern.candidate_star_v_maps[star_idx].size(); i++) {
            // map candidate_star_v_maps[idx][i] to star_vertexes[idx]
            for (uint32_t j = 0; j < pattern.candidate_star_v_maps[star_idx][i].size(); j++) {
                v_map.emplace(pattern.candidate_star_v_maps[star_idx][i][j], star_vertexes[star_idx][j]);
            }

            // DFS_judge_further
            DFS_judge_remaining_edges(idx + 1, graph, pattern, cycle_vertexes, star_vertexes, v_map, found_instances, lost_edge_cases);

            // restore
            for (uint32_t j = 0; j < pattern.candidate_star_v_maps[star_idx][i].size(); j++) {
                v_map.erase(pattern.candidate_star_v_maps[star_idx][i][j]);
            }
        }
    }
}


vector<uint64_t> judge_remaining_edges(Graph &graph, vector<Pattern> &patterns,
                               vector <vector<uint32_t>> &cycle_vertexes, vector <vector<uint32_t>> &star_vertexes, uint32_t *lost_edge_cases) {
    /*
    for (uint32_t i = 0; i < cycle_vertexes.size(); i++) {
        for (uint32_t j = 0; j < cycle_vertexes[i].size(); j++)
            cout << cycle_vertexes[i][j] << " ";
        cout << endl;
    }

    for (uint32_t i = 0; i < star_vertexes.size(); i++) {
        for (uint32_t j = 0; j < star_vertexes[i].size(); j++)
            cout << star_vertexes[i][j] << " ";
        cout << endl;
    } */

    vector<uint64_t> found_instances(patterns.size(), 0);
    for (uint32_t k = 0; k < patterns.size(); k++) {
        map <uint32_t, uint32_t> v_map;
        DFS_judge_remaining_edges(0, graph, patterns[k], cycle_vertexes, star_vertexes, v_map, found_instances[k], lost_edge_cases);
    }
    
    return found_instances;
}

/* Major helper function to iterate a sampling tree
*/
void tree_DFS_helper(SamplingTreeNode *root, Graph &graph, vector<Pattern> &patterns,
                     vector<SamplingTreeNode *> &path, map <uint32_t, uint32_t> &v_map, vector<double> &total_prob_inverse, uint64_t & total_sampling_times,
                     uint32_t *lost_edge_cases) {
    // corner cases
    if (root->children.size() == 0) {
        path.push_back(root);
        v_map.clear();

        // task 1: print path
        assert(path.size() == (patterns[0].cycle_num + patterns[0].star_num) * 2); // since all patterns have same star_num and cycle_num and thus node_num

        if (FLAG_print) {
            cout << "found a path: ";
        }

        set <uint32_t> vertex_set;
        bool test_succeed = true;
        double prob_inverse = 1.0;
        total_sampling_times += 1;

        // task 2: judge cycles
        vector <vector<uint32_t>> cycle_vertexes;
        for (uint32_t i = 0; i < patterns[0].cycle_num; i++) {
            set <uint32_t> cycle_vertex_set;
            SamplingTreeNode *cur_root = path[i * 2];
            SamplingTreeNode *cur_leave = path[i * 2 + 1];

            if (FLAG_print) {
                cur_root->print();
                cout << ", ";
                cur_leave->print();
                cout << ", ";
            }

            // handle cycle root
            prob_inverse *= cur_root->probability;

            // map this node's vertexes and cycles[i]'s first several nodes
            vector <uint32_t> cur_cycle;
            for (uint32_t j = 0; j < cur_root->vertexes.size(); j++) {
                v_map.emplace(patterns[0].cycles[i][j], cur_root->vertexes[j]);
                vertex_set.insert(cur_root->vertexes[j]);
                cur_cycle.push_back(cur_root->vertexes[j]);
                cycle_vertex_set.insert(cur_root->vertexes[j]);
            }

            // handle cycle leave
            prob_inverse *= cur_leave->probability;

            v_map.emplace(patterns[0].cycles[i].back(), cur_leave->vertexes[0]);
            vertex_set.insert(cur_leave->vertexes[0]);
            cur_cycle.push_back(cur_leave->vertexes[0]);
            cycle_vertexes.push_back(cur_cycle);
            cycle_vertex_set.insert(cur_leave->vertexes[0]);

            // judge cycle construct
            // test overlapped vertex
            if (cycle_vertex_set.size() != patterns[0].cycles[i].size()) {
                test_succeed = false;
                if (FLAG_print)
                    cout << "cycle " << i << " has overlapped vertexes" << endl;
            }

            // test copy
            if (!test_circle_copy(graph, cur_root->vertexes, cur_leave->vertexes[0])) {
                test_succeed = false;
                if (FLAG_print)
                    cout << "cycle " << i << " does not pass copy test" << endl;
            }
        }

        // task 3: judge stars
        vector <vector<uint32_t>> star_vertexes;
        if (test_succeed != false) {
            for (uint32_t i = 0; i < patterns[0].star_num; i++) {
                SamplingTreeNode *cur_root = path[(patterns[0].cycle_num + i) * 2];
                SamplingTreeNode *cur_leave = path[(patterns[0].cycle_num + i) * 2 + 1];

                if (FLAG_print) {
                    cur_root->print();
                    cout << ", ";
                    cur_leave->print();
                    cout << ", ";
                }

                // handle star root
                prob_inverse *= cur_root->probability;
                v_map.emplace(patterns[0].stars[i].front(), cur_root->vertexes[0]);
                vector <uint32_t> cur_star;
                cur_star.push_back(cur_root->vertexes[0]);
                vertex_set.insert(cur_root->vertexes[0]);

                // handle leaves and and test copy
                prob_inverse *= cur_leave->probability;
                uint32_t petal_num = cur_leave->vertexes.size();
                if (petal_num == 0) {
                    test_succeed = false;
                    if (FLAG_print)
                        cout << "star " << i << " does not have enough petals" << endl;
                } /* else if (petal_num == 1 && !graph.vertex_before(cur_root->vertexes[0], cur_leave->vertexes[0])) {
                    test_succeed = false;
                    if (FLAG_print)
                        cout << "star " << i << " does not form a copy" << endl;
                } */
                else {
                    for (uint32_t j = 0; j < cur_leave->vertexes.size(); j++) {
                        v_map.emplace(patterns[0].stars[i][j + 1], cur_leave->vertexes[j]);
                        cur_star.push_back(cur_leave->vertexes[j]);
                        vertex_set.insert(cur_leave->vertexes[j]);
                    }
                    star_vertexes.push_back(cur_star);
                }
            }
            if (FLAG_print)
                cout << endl;
        }

        // print v_map
        // assert(v_map.size() == pattern.vertex_count);

        // task 4: judge other edges
        vector<uint64_t> num_matches(patterns.size(), 0);
        if (test_succeed == false) {
            if (FLAG_print)
                cout << "previous cycle/star tests already failed" << endl;
        } else if (vertex_set.size() != patterns[0].vertex_count) {
            if (FLAG_print)
                cout << "this path has overlapped sampled vertexes" << endl;
            test_succeed = false;
        } else {
            // judge other edges
            //if (!judge_remaining_edges(graph, pattern, cycle_vertexes, star_vertexes)) {
            //    test_succeed = false;
            //}
            num_matches = judge_remaining_edges(graph, patterns, cycle_vertexes, star_vertexes, lost_edge_cases);
        }

        // task 4: if pass test, count
        if (test_succeed) {
            if (FLAG_print) {
                cout << "the path passed the test!" << endl;
            }
            for (uint32_t k = 0; k < patterns.size(); k++)
                total_prob_inverse[k] += prob_inverse * num_matches[k];
        }
    

        path.pop_back();
        return;
    }

    // DFS to iterate children
    path.push_back(root);
    for (uint32_t i = 0; i < root->children.size(); i++) {
        tree_DFS_helper(root->children[i], graph, patterns, path, v_map, total_prob_inverse, total_sampling_times, lost_edge_cases);
    }
    path.pop_back();

    return;
}

/* Count a sampling tree
*/
void sampling_tree_count(SamplingTreeNode *root, Graph &graph, vector<Pattern> &patterns, vector<double> &total_prob_inverse, uint64_t & total_sampling_times, uint32_t *lost_edge_cases) {
    if (FLAG_print)
        cout << "====== start to count this tree" << endl;

    // DFS iterate from root to a leave
    vector < SamplingTreeNode * > path;
    map <uint32_t, uint32_t> v_map;

    tree_DFS_helper(root, graph, patterns, path, v_map, total_prob_inverse, total_sampling_times, lost_edge_cases);

    return;
}

void free_sampling_tree(SamplingTreeNode *root)
{
    if (root->children.size() == 0)
    {
        delete root;
        root = nullptr;
        return ;
    }

    for (uint32_t i = 0; i < root->children.size(); i++)
    {
        free_sampling_tree(root->children[i]);
    }
    
    delete root;
    root = nullptr;
}

/* Estimating arbitrary patterns
* from graph, with sampling_times
*/
tuple<vector<double>, uint64_t> estimate_pattern(uint32_t thread_id, Graph &graph, vector<Pattern> &patterns, uint64_t sampling_times, uint32_t *lost_edge_cases, bool use_cache, default_random_engine &rand_generator) {
    //cout << "= Start to estimate patterns, " << sampling_times << " times" << endl;
    vector<double> total_prob_inverse(patterns.size(), 0.0);

    for (uint64_t i = 0; i < sampling_times; i++) {
        // Step 1: construct the trees
        SamplingTreeNode *root = construct_sampling_tree(thread_id, graph, patterns[0], use_cache, rand_generator);

        // Step 2: iterate every path from root to leave
        vector<double> this_prob_inverse(patterns.size(), 0.0);
        uint64_t this_num_leaves = 0;
        sampling_tree_count(root, graph, patterns, this_prob_inverse, this_num_leaves, lost_edge_cases);
        
        // Step 3: free sampling tree
        free_sampling_tree(root);


        for (uint32_t k = 0; k < patterns.size();k++)
        {   
            total_prob_inverse[k] += this_prob_inverse[k] / this_num_leaves;
        }

    }

    return make_tuple(total_prob_inverse, sampling_times);
}

void init_cache_index() {
    for (uint32_t i = 0; i < MAX_THREAD_NUM; i++) {
        for (auto it = cycle_cache_index[i].begin(); it != cycle_cache_index[i].end(); it++)
            it->second = 0;
        for (auto it = star_cache_index[i].begin(); it != star_cache_index[i].end(); it++)
            it->second = 0;
    }
}

void append_next_cache() {
    for (uint32_t i = 0; i < MAX_THREAD_NUM; i++) {   
        for (auto it = cycle_cache_next[i].begin(); it != cycle_cache_next[i].end(); it++) {
            if (cycle_cache[i].find(it->first) == cycle_cache[i].end()) {
                cycle_cache[i][it->first] = it->second;
                continue;
            }
            vector<struct SamplingTreeNode *> & vec = cycle_cache[i].at(it->first);
            vec.insert(vec.end(), it->second.begin(), it->second.end());
        }
        for (auto it = star_cache_next[i].begin(); it != star_cache_next[i].end(); it++) {
            if (star_cache[i].find(it->first) == star_cache[i].end()) {
                star_cache[i][it->first] = it->second;
                continue;
            }
            vector<struct SamplingTreeNode *> & vec = star_cache[i].at(it->first);
            vec.insert(vec.end(), it->second.begin(), it->second.end());
        }
    }
}

#endif
