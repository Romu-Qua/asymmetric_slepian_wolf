//
// Created by bob on 8/3/22.
//

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <random>
#include <tuple>
#include "simulation_utils.h"

#ifndef INFORMATION_THEORY_ENCODING_DECODING_H
#define INFORMATION_THEORY_ENCODING_DECODING_H

#endif //INFORMATION_THEORY_ENCODING_DECODING_H

using namespace std;

tuple<bool, vector<bool>> decode_at_current_rate(const vector<double> &llrs,
                                                 const vector<bool> &syndrome,
                                                 const int n_cols,
                                                 const int n_rows,
                                                 const vector<vector<int>> &pos_varn,
                                                 const vector<vector<int>> &pos_checkn,
                                                 const vector<uint32_t> column_pointers,
                                                 const vector<uint16_t> row_index,
                                                 const std::size_t max_num_iter,
                                                 const double vsat);

template<typename T>
static void saturate(vector<vector<T>> &mv, const double vsat);


void hard_decision(vector<bool> &out, const vector<double> &llrs, const vector<vector<double>> &msg_c);
void var_node_update(vector<vector<double>> &msg_v,
                     const vector<vector<double>> &msg_c,
                     const vector<double> &llrs,
                     const vector<vector<int>> &pos_varn,
                     const vector<vector<int>> &pos_checkn,
                     const int n_cols);

void check_node_update(vector<vector<double>> &msg_c,
                       const vector<vector<double>> &msg_v,
                       const vector<bool> &syndrome,
                       const vector<vector<int>> &pos_varn,
                       const vector<vector<int>> &pos_checkn,
                       const int n_cols,
                       const int n_rows);

tuple<vector<vector<int>>, vector<vector<int>>> calculate_vn_cv(int n_cols,
                                                                int n_rows,
                                                                vector<uint32_t> column_pointers,
                                                                vector<uint16_t> row_index);

vector<double> bsc_llr(vector<bool> &y, double p);

vector<bool>  encode(vector<bool> &in, vector<uint32_t> column_pointers,
                     vector<uint16_t> row_index, uint16_t n_rows);
