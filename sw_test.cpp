//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <cstdint>
#include <iostream>
#include "simulation_utils.h"
#include "encoding_decoding.h"
#include "npy.hpp"

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

This is a test to figure out how to use c++,
it follows the code from here https://github.com/XQP-Munich/LDPC4QKD
to a great extend



 This is the main file of this project, it contains the read in of the
 matrices, the encoding and decoding of the data and the calculation of the
 FER, as well as the storage of the results.
 If you have different codes or want to change paramters, do it here.

*/





using namespace std;

#ifdef DEBUG_MESSAGES_ENABLED

#include <iostream>

#define DEBUG_MESSAGE(msg) do {std::cerr << msg << std::endl;} while (false)

#else

#define DEBUG_MESSAGE(msg)

#endif /* ifdef DEBUG_MESSAGES_ENABLED */


// example of how to express your example code directly in this file instead of using
// numpy
///    H =  [1 0 1 0 1 0 1
    ///			0 1 1 0 0 1 1
    ///			0 0 0 1 1 1 1]

/// writing the matrix in csc format
//vector<uint32_t> column_pointers{0, 1, 2, 4, 5, 7, 9, 12};
//vector<uint16_t> row_index{0, 1, 0, 1, 2, 0, 2, 1, 2, 0, 1, 2};


// all the parameters to use different codes, and to set the location of where to save the results
uint16_t n_rows = 212;
uint16_t n_cols = 1908;
double sweep_min = 0.01088889;
double sweep_max = 0.01444445;
int sweep_steps = 5;
const char * path = {"codes/1908_212_4_colmn_pointers.npy"};
const char * path2 = {"codes/1908_212_4_row_index.npy"};
string path_fer("results/fer_detail_1908_212_4_big_error");
string path_p("results/p_detail_1908_212_4_big_error");
int number_of_samples = 100;

// templates in relation to numpy arrays
template <typename Scalar>
struct npy_data {
    std::vector<unsigned long> shape;
    bool fortran_order;
    std::vector<Scalar> data;
};

template <typename Scalar>
npy_data<Scalar> test_load(const char * path) {
    npy_data<Scalar> d;
    npy::LoadArrayFromNumpy(path, d.shape, d.fortran_order, d.data);
    return d;
}


/** Simulates a BSC channel with a given error probability
 *
 *  @param data_size block size of the code/ message length
 *  @param p crossover probability of the BSC
 *  @param column_pointers column pointers of the code
 *  @param row_index row indeces of the code
 *  @return True if decoded succesfully before max iterations run out, false otherwise
 */
bool simulate_bsc(const int data_size,
                  const double p,
                  const vector<uint32_t> column_pointers,
                  const vector<uint16_t> row_index){

    // generate random input
    vector<bool> input = random_input(data_size);
    vector<bool> checksum;
    checksum = encode(input, column_pointers, row_index, n_rows);

    //cout << "encoded \n";

    vector<bool> y = bit_flip_channel(input, p);

    // calculating the log-likehood ratios
    vector<double> llr_init = bsc_llr(y, p);

    auto [pos_varn, pos_checkn] = calculate_vn_cv(n_cols, n_rows, column_pointers, row_index);

    // decoding, max iterations is at 50 right now
    auto  [success, x_prime] = decode_at_current_rate(llr_init, checksum, n_cols, n_rows, pos_varn, pos_checkn,
                                                      column_pointers, row_index, 50, 100);

    //cout << "decoded \n" << "with success: " << success << endl;
    //cout << "vectors are equal: " << (x_prime == input) << endl;

    bool true_success = x_prime == input;
    return true_success;
}

/** main function starting the simulation and saving the results
 */
int main() {

    vector<double> p_vec = linspace(sweep_min, sweep_max, sweep_steps);
    vector<double> fers = vector<double>(p_vec.size());
    int data_size = n_cols;
    int number_of_success = 0;
    double p;

    // loading the code from the numpy arrays
    auto d = test_load<unsigned int>(path);
    auto d2 = test_load<uint16_t>(path2);
    vector<uint32_t>column_pointers = d.data;
    vector<uint16_t>row_index = d2.data;

    // looping over all samples
    for (size_t i = 0; i < p_vec.size(); ++i) {
        number_of_success = 0;
        p = p_vec[i];
        for (int i = 0; i < number_of_samples; i++) {
            cout << "running sample: " << i << endl;
            if (simulate_bsc(data_size, p, column_pointers, row_index)) {
                number_of_success++;
            }
        }
        cout << "number of success: " << number_of_success << endl;
        fers[i] = (number_of_samples - (double) number_of_success) / number_of_samples;
        cout << "current frame error rate: " << fers[i] << "for ber " << p << endl;
    }

    //cout << "number of success: " << number_of_success << endl;
    //cout << "frame error rate: " << (number_of_samples-(double)number_of_success) / number_of_samples << endl;
    print(fers);

    // saving the results
    ofstream file_fer(path_fer, ios::out | ofstream::binary);
    copy(fers.begin(), fers.end(), ostream_iterator<double>(file_fer, "\n"));
    file_fer.close();


    ofstream file_p(path_p, ios::out | ofstream::binary);
    copy(p_vec.begin(), p_vec.end(), ostream_iterator<double>(file_p, "\n"));
    file_p.close();

    return 0;
}


