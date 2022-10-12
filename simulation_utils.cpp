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


/* 
    * @brief: applies a bit flip with probability p
    * @param: input data, probability of bit flip
    * @return: noisy data
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

/*
    * @brief: applies a bit flip with probability p
    * @param: input data, probability of bit flip
    * @return: noisy data
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


void print(vector<bool> const &input)
{
    for (auto it = input.cbegin(); it != input.cend(); it++){
        cout << *it << ' ';
    }
    cout << endl;
}


void print(vector<double> const &input)
{
    for (auto it = input.cbegin(); it != input.cend(); it++){
        cout << *it << ' ';
    }
}


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



int number_of_diff_vector_elements(vector<bool> const &input1, vector<bool> const &input2){
    int count = 0;
    for (int i = 0; i < input1.size(); i++){
        if (input1[i] != input2[i]){
            count += 1;
        }
    }
    return count;
}


vector<bool> random_input(const int size){
    static mt19937 generator;
    static uniform_int_distribution<int> distribution(0,1);
    vector<bool> vec(size, false);

    for(auto val : vec) val = distribution(generator);

    return vec;
}

void dummy(){
    std::cout << "idiot" << std::endl;
}
