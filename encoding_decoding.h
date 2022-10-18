/**
 * Copyright (c) 2022 Ronny Mueller ronny.r_mueller@web.de

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is furnished
to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


THIs code follows at takes from here https://github.com/XQP-Munich/LDPC4QKD
to a great extend

header file for the actual decoder implementation
*/


//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
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



/**
 * @brief Tries to decode the given codeword using the given parity check matrix
 * @param llrs inital log-likelihood ratios
 * @param syndrome The checksum/syndrome of the codeword
 * @param n_cols number of columns of H
 * @param n_rows number of rows of H
 * @param pos_varn List of lists of positions of variable nodes
 * @param pos_checkn List of lists of positions of check nodes
 * @param column_pointers column pointers of H in CSC
 * @param row_index row indices of H in CSC
 * @param max_num_iter max number of decoding iterations
 * @param vsat cut-off value for messages
 * @return
 */
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



/**
 * @brief caps the values of a vector in both positive and negative direction
 * @tparam T
 * @param mv The vector to cap
 * @param vsat the limit to cap, +/-
 */
template<typename T>
static void saturate(vector<vector<T>> &mv, const double vsat);

/**
 * @brief sums up all messages to calculate if the llr is negative or positive, returns
 * the current most likely bit
 * @param out the vector to be fille dwith the bits
 * @param llrs the initial likelihoods
 * @param msg_c the current messages
 */
void hard_decision(vector<bool> &out, const vector<double> &llrs, const vector<vector<double>> &msg_c);



/**
 * @brief performs the variable node update step
 * @param msg_v the array containing the variable messages
 * @param msg_c the array containing the check messages
 * @param llrs intial log likelihood ratios
 * @param pos_varn List of lists of positions of variable nodes
 * @param pos_checkn List of lists of positions of check nodes
 * @param n_cols number of columns of H
 */
void var_node_update(vector<vector<double>> &msg_v,
                     const vector<vector<double>> &msg_c,
                     const vector<double> &llrs,
                     const vector<vector<int>> &pos_varn,
                     const vector<vector<int>> &pos_checkn,
                     const int n_cols);


/**
 * @brief performs the check node update step
 * @param msg_c the array containing the check messages
 * @param msg_v the array containing the variable messages
 * @param syndrome the syndrome of the codeword
 * @param pos_varn List of lists of positions of variable nodes
 * @param pos_checkn List of lists of positions of check nodes
 * @param n_cols number of columns of H
 * @param n_rows number of rows of H
 */
void check_node_update(vector<vector<double>> &msg_c,
                       const vector<vector<double>> &msg_v,
                       const vector<bool> &syndrome,
                       const vector<vector<int>> &pos_varn,
                       const vector<vector<int>> &pos_checkn,
                       const int n_cols,
                       const int n_rows);


/**
 * @brief  represent the matrix in terms of list of lists from row and column view
 * @param n_cols number of columns of H
 * @param n_rows number of rows of H
 * @param column_pointers column pointers of H in CSC
 * @param row_index row indices of H in CSC
 * @return List of list of positions of variable nodes and List of list of positions of check nodes
 */
tuple<vector<vector<int>>, vector<vector<int>>> calculate_vn_cv(int n_cols,
                                                                int n_rows,
                                                                vector<uint32_t> column_pointers,
                                                                vector<uint16_t> row_index);

/**
 * @brief calculates the initial log likelihood ratios
 * @param y the received message
 * @param p the crossover probability
 * @return log-likelihood rations
 */
vector<double> bsc_llr(vector<bool> &y, double p);


/**
 * @brief calculates the syndrome of the codeword (matrix multiplication of H and y in CSC)
 * @param in the input vector to multply with H
 * @param column_pointers column pointers of H in CSC
 * @param row_index row indices of H in CSC
 * @param n_rows number of rows of H
 * @return the matrix product (the syndrome)
 */
vector<bool>  encode(vector<bool> &in, vector<uint32_t> column_pointers,
                     vector<uint16_t> row_index, uint16_t n_rows);
