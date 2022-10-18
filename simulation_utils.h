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



#include <vector>


/**
 * @brief This function applies a bit flip to a given vector with probability p
 *
 * @param in The vector to which the bit flip is applied
 * @param p The probability of the bit flip
 *
 * @return The vector with the bit flip applied
 */
std::vector<bool> bit_flip_channel(std::vector<bool> in, double p);



/**
 * @brief Applies a specific number of bit flip error to vector
 * @param in The vector to which the bit flip is applied
 * @param number_of_errors  The number of bit flips to apply
 * @return The vector with the bit flip applied
 */
std::vector<bool> bit_flip_channel_det(std::vector<bool> &in, int number_of_errors);



/**
 * @brief prints a vector of bools
 *
 * @param input The vector to print
 */
void print(std::vector<bool> const &input);



/**
 * @brief prints a vector of doubles
 *
 * @param input The vector to print
 */
void print(std::vector<double> const &input);



/**
 * @brief creates a vector with linear steps between start and end
 * @param min start of the vector
 * @param max end of the vector
 * @param n number of steps
 * @return vector of linear steps
 */
std::vector<double> linspace(double min, double max, int n);


/**
 * @brief calculates the number of differences between two vectors
 * @param input1 vector 1
 * @param input2 vector 2
 * @return number of differences
 */
int number_of_diff_vector_elements(std::vector<bool> const &input1, std::vector<bool> const &input2);

/**
 * @brief creates a random vector of bools
 * @param size size of the vector
 * @return random vector of bools
 */
std::vector<bool> random_input(const int size);




