#include <vector>
#include <cmath>
#include <cstdint>
#include <iostream>
#include "simulation_utils.h"
#include "encoding_decoding.h"
#include "npy.hpp"

/// This is a test to figure out how to use c++,
/// it follows the code from here https://github.com/XQP-Munich/LDPC4QKD

using namespace std;
//using namespace simulation_utils;



#ifdef DEBUG_MESSAGES_ENABLED

#include <iostream>

#define DEBUG_MESSAGE(msg) do {std::cerr << msg << std::endl;} while (false)

#else

#define DEBUG_MESSAGE(msg)

#endif /* ifdef DEBUG_MESSAGES_ENABLED */

///    H =  [1 0 1 0 1 0 1
    ///			0 1 1 0 0 1 1
    ///			0 0 0 1 1 1 1]

/// writing the matrix in csc format
//vector<uint32_t> column_pointers{0, 1, 2, 4, 5, 7, 9, 12};
//vector<uint16_t> row_index{0, 1, 0, 1, 2, 0, 2, 1, 2, 0, 1, 2};
uint16_t n_rows = 738;
uint16_t n_cols = 4095;

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

    // decoding
    auto  [success, x_prime] = decode_at_current_rate(llr_init, checksum, n_cols, n_rows, pos_varn, pos_checkn,
                                                      column_pointers, row_index, 50, 100);

    //cout << "decoded \n" << "with success: " << success << endl;
    //cout << "vectors are equal: " << (x_prime == input) << endl;

    bool true_success = x_prime == input;
    return true_success;
}


int main() {

    vector<double> p_vec = linspace(0.005, .03, 10);
    vector<double> fers = vector<double>(p_vec.size());
    int data_size = 4095;
    double p = 0.04;
    int number_of_samples = 1000;
    int number_of_success = 0;

    // I haven't quite figured out how to use relative paths here
    const char * path = {"/home/bob/Documents/information_theory/4095_738_102_colmn_pointers.npy"};
    const char * path2 = {"/home/bob/Documents/information_theory/4095_738_102_row_index.npy"};



    //const char * path = {"/home/bob/Documents/information_theory/4095_737_101_colmn_pointers.npy"};
    //const char * path2 = {"/home/bob/Documents/information_theory/4095_737_101_row_index.npy"};
    auto d = test_load<unsigned int>(path);
    auto d2 = test_load<uint16_t>(path2);
    vector<uint32_t>column_pointers = d.data;
    vector<uint16_t>row_index = d2.data;

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

    string path_fer("/home/bob/Documents/information_theory/fer_detail_4095_738_102");
    ofstream file_fer(path_fer, ios::out | ofstream::binary);
    copy(fers.begin(), fers.end(), ostream_iterator<double>(file_fer, "\n"));
    file_fer.close();

    string path_p("/home/bob/Documents/information_theory/p_detail_4095_738_102");
    ofstream file_p(path_p, ios::out | ofstream::binary);
    copy(p_vec.begin(), p_vec.end(), ostream_iterator<double>(file_p, "\n"));
    file_p.close();

    return 0;
}


