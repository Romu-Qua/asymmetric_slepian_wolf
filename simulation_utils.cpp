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

This file just contains a bunch of utility functions for the simulation
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

using namespace std;

#ifdef DEBUG_MESSAGES_ENABLED

#include <iostream>

#define DEBUG_MESSAGE(msg) do {std::cerr << msg << std::endl;} while (false)

#else

#define DEBUG_MESSAGE(msg)

#endif /* ifdef DEBUG_MESSAGES_ENABLED */


/**
 * @brief This function applies a bit flip to a given vector with probability p
 *
 * @param in The vector to which the bit flip is applied
 * @param p The probability of the bit flip
 *
 * @return The vector with the bit flip applied
 */
vector<bool> bit_flip_channel(vector<bool> in, double p) {
    srand(time(NULL));
    vector<bool> out = in;

    random_device rd;
    mt19937 gen(rd());
    bernoulli_distribution d(p);
    for (size_t i = 0; i<in.size(); i++) {
        if (d(gen)) {
            out[i] = !out[i];
        }
    }
    return out;
}

/**
 * @brief Applies a specific number of bit flip error to vector
 * @param in The vector to which the bit flip is applied
 * @param number_of_errors  The number of bit flips to apply
 * @return The vector with the bit flip applied
 */
vector<bool> bit_flip_channel_det(vector<bool> &in, int number_of_errors) {
    vector<bool> out = in;
    vector<int> used_idx = vector<int>(number_of_errors, -1);
    srand(time(NULL));
    // go over the number of required flips
    for (size_t i = 0; i<number_of_errors; i++){

        // random index
        size_t index = rand() % in.size();
        // if this index flipped already, resample until it is not flipped
        while (find(used_idx.begin(), used_idx.end(), index) != used_idx.end()){
            index = rand() % in.size();
        }
        // append flipped to list of flipped
        used_idx[i] = index;
        out[index] = !out[index];
    }
    return out;
}

/**
 * @brief prints a vector of bools
 *
 * @param input The vector to print
 */
void print(vector<bool> const &input)
{
    for (auto it = input.cbegin(); it != input.cend(); it++){
        cout << *it << ' ';
    }
    cout << endl;
}

/**
 * @brief prints a vector of doubles
 *
 * @param input The vector to print
 */
void print(vector<double> const &input)
{
    for (auto it = input.cbegin(); it != input.cend(); it++){
        cout << *it << ' ';
    }
}

/**
 * @brief creates a vector with linear steps between start and end
 * @param min start of the vector
 * @param max end of the vector
 * @param n number of steps
 * @return vector of linear steps
 */
vector<double> linspace(double min, double max, int n)
{
    vector<double> result;
    // vector iterator
    int iterator = 0;

    for (int i = 0; i <= n-2; i++)
    {
        double temp = min + i*(max-min)/(floor((double)n) - 1);
        result.insert(result.begin() + iterator, temp);
        iterator += 1;
    }

    //iterator += 1;

    result.insert(result.begin() + iterator, max);
    return result;
}


/**
 * @brief calculates the number of differences between two vectors
 * @param input1 vector 1
 * @param input2 vector 2
 * @return number of differences
 */
int number_of_diff_vector_elements(vector<bool> const &input1, vector<bool> const &input2){
    int count = 0;
    for (int i = 0; i < input1.size(); i++){
        if (input1[i] != input2[i]){
            count += 1;
        }
    }
    return count;
}

/**
 * @brief creates a random vector of bools
 * @param size size of the vector
 * @return random vector of bools
 */
vector<bool> random_input(const int size){
    static mt19937 generator;
    static uniform_int_distribution<int> distribution(0,1);
    vector<bool> vec(size, false);

    for(auto val : vec) val = distribution(generator);

    return vec;
}

