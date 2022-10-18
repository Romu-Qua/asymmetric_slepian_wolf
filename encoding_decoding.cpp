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
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>
#include <tuple>
#include "encoding_decoding.h"

using namespace std;




#ifdef DEBUG_MESSAGES_ENABLED

#include <iostream>

#define DEBUG_MESSAGE(msg) do {std::cerr << msg << std::endl;} while (false)

#else

#define DEBUG_MESSAGE(msg)

#endif /* ifdef DEBUG_MESSAGES_ENABLED */



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
                                                 const std::size_t max_num_iter = 50,
                                                 const double vsat = 100) {
    // check inputs.
    if (llrs.size() != n_cols) {
        throw runtime_error("input doesn't match H.");
    }

    if (syndrome.size() != n_rows) {
        throw runtime_error(
                "checksum doesn't match number of rows in H");
    }

    vector<bool> out = vector<bool>(llrs.size(), false);

    vector<vector<double>> msg_v(n_rows);  // messages from variable nodes to check nodes
    vector<vector<double>> msg_c(n_cols);  // messages from check nodes to variable nodes

    // initialize msg_v
    for (size_t i{}; i < msg_v.size(); ++i) {
        auto &curr_mv = msg_v[i];
        curr_mv.resize(pos_varn[i].size());
        for (size_t j{}; j < msg_v[i].size(); ++j) {
            curr_mv[j] = llrs[pos_varn[i][j]];
        }
    }

    // initialize msg_c
    for (size_t i{}; i < msg_c.size(); ++i) {
        msg_c[i].resize(pos_checkn[i].size());
    }

    for (size_t it_unused{}; it_unused < max_num_iter; ++it_unused) {
        check_node_update(msg_c, msg_v, syndrome, pos_varn, pos_checkn, n_cols, n_rows);
        saturate(msg_c, vsat);

        var_node_update(msg_v, msg_c, llrs, pos_varn, pos_checkn, n_cols);
        saturate(msg_v, vsat);

        // hard decision
        hard_decision(out, llrs, msg_c);

        // terminate decoding if codeword matches syndrome
        vector<bool> decision_syndrome(syndrome.size());
        decision_syndrome = encode(out, column_pointers, row_index, n_rows);
        if (decision_syndrome == syndrome) {
            return {true, out};
        }

        // check for diverging decoder
        for (const auto &m : msg_v) {
            for (const auto &v : m) {
                if (std::isnan(v)) {
                    return {false, out};
                }
            }
        }
    }

    return {false, out};  // Decoding was not successful.
}


/**
 * @brief caps the values of a vector in both positive and negative direction
 * @tparam T
 * @param mv The vector to cap
 * @param vsat the limit to cap, +/-
 */
template<typename T>
static void saturate(vector<vector<T>> &mv,
                     const double vsat) {
    for (auto &v : mv) {
        for (auto &a : v) {
            if (a > vsat) { a = vsat; }
            else if (a < -vsat) { a = -vsat; }
        }
    }
}


/**
 * @brief calculates the syndrome of the codeword (matrix multiplication of H and y in CSC)
 * @param in the input vector to multply with H
 * @param column_pointers column pointers of H in CSC
 * @param row_index row indices of H in CSC
 * @param n_rows number of rows of H
 * @return the matrix product (the syndrome)
 */
vector<bool>  encode(vector<bool> &in, vector<uint32_t> column_pointers,
                     vector<uint16_t> row_index, uint16_t n_rows) {


    // initialize the output vector to all 0/false
    vector<bool> out = vector<bool>(n_rows, false);

    // iterate over the input vector
    for (size_t col = 0; col<in.size(); col++){
        for (size_t j = column_pointers[col]; j<column_pointers[col+1]; j++){
            out[row_index[j]] = out[row_index[j]] != in[col];

        }
    }
    return out;
}


/**
 * @brief calculates the initial log likelihood ratios
 * @param y the received message
 * @param p the crossover probability
 * @return log-likelihood rations
 */
vector<double> bsc_llr(vector<bool> &y, double p) {
    vector<double> llr = vector<double>(y.size(), 0);
    for (size_t i = 0; i<y.size(); i++){
        if (y[i]){
            llr[i] = log(p/(1-p));
        }
        else{
            llr[i] = log((1-p)/p);
        }
    }
    //cout << y.size() << endl;
    return llr;
}


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
                                                                vector<uint16_t> row_index) {

    vector<vector<int>> pos_varn = vector<vector<int>>(n_rows, vector<int>{});
    vector<vector<int>> pos_checkn = vector<vector<int>>(n_cols, vector<int>{});

    for (int col = 0; col < n_cols; col++) {
        for (size_t j = column_pointers[col]; j < column_pointers[col + 1u]; j++) {
            pos_varn[row_index[j]].push_back(col); // I am pretty sure this is terrible, new alloctiaon each iter
        }
    }

    pos_checkn.assign(n_cols, vector<int>{});

    for (int i{}; i < pos_varn.size(); ++i) {
        for (auto &vn : pos_varn[i]) {
            pos_checkn[vn].push_back(i);
        }
    }
    return {pos_varn, pos_checkn};
}


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
                       const int n_rows) {
    double msg_part{};
    vector<size_t> mc_position(n_cols);

    for (size_t m{}; m < n_rows; ++m) {
        // product of incoming messages
        double mc_prod = 1 - 2 * static_cast<double>(syndrome[m]);
        // Note: pos_varn[m].size() = check_node_degrees[m]
        const auto curr_check_node_degree = pos_varn[m].size();
        for (std::size_t k{}; k < curr_check_node_degree; ++k) {
            mc_prod *= ::tanh(0.5 * msg_v[m][k]);
        }

        for (std::size_t k{}; k < curr_check_node_degree; ++k) {
            // computing message from
            if (msg_v[m][k] == 0.) {
                msg_part = 1;
                for (std::size_t non_k{}; non_k < curr_check_node_degree; ++non_k) {
                    if (non_k != k) {
                        msg_part *= ::tanh(0.5 * msg_v[m][k]);
                    }
                }
            } else {
                msg_part = mc_prod / ::tanh(0.5 * msg_v[m][k]);
            }

            auto msg_final = ::log((1 + msg_part) / (1 - msg_part));

            // place the message at the correct position in the output array
            const int curr_pos_varn = pos_varn[m][k];
            msg_c[curr_pos_varn][mc_position[curr_pos_varn]] = msg_final;
            mc_position[curr_pos_varn]++;
        }
    }
}


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
                     const int n_cols){
    vector<size_t> mv_position(n_cols);

    for (size_t m{}; m < llrs.size(); ++m) {
        const double mv_sum = accumulate(msg_c[m].begin(), msg_c[m].end(), llrs[m]);

        // Note: pos_checkn[m].size() = var_node_degs[m]
        for (size_t k{}; k < pos_checkn[m].size(); ++k) {
            const double msg = mv_sum - msg_c[m][k];

            // place the message at the correct position in the output array
            const int curr_pos_cn = pos_checkn[m][k];
            msg_v[curr_pos_cn][mv_position[curr_pos_cn]] = msg;
            mv_position[curr_pos_cn]++;
        }
    }
}


/**
 * @brief sums up all messages to calculate if the llr is negative or positive, returns
 * the current most likely bit
 * @param out the vector to be fille dwith the bits
 * @param llrs the initial likelihoods
 * @param msg_c the current messages
 */
void hard_decision(
        vector<bool> &out,
        const vector<double> &llrs,
        const vector<vector<double>> &msg_c){
    fill(out.begin(), out.end(), 0);
    for (size_t j{}; j < llrs.size(); ++j) {
        const double curr_sum = accumulate(msg_c[j].begin(), msg_c[j].end(), llrs[j]);
        if (curr_sum < 0) {
            out[j] = 1;
        }
    }
}





