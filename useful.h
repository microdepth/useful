#ifndef USEFUL_USEFUL_H
#define USEFUL_USEFUL_H

#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <math.h>
using namespace std;

// -----------------------------------------------------------------------------
// 1. CLASSES
// -----------------------------------------------------------------------------

// Class to handle teeing (duplicating) the output stream to two different streams.
class TeeBuf : public streambuf {
public:
    TeeBuf(streambuf* sb1, streambuf* sb2) : buf1(sb1), buf2(sb2) {};

    virtual int overflow(int c) {
        if (c != EOF) {
            buf1->sputc(c);
            buf2->sputc(c);
        }
        return c;
    }

    virtual int sync() {
        buf1->pubsync();
        buf2->pubsync();
        return 0;
    }
private:
    streambuf* buf1;
    streambuf* buf2;
};

// -----------------------------------------------------------------------------
// 2. MATHEMATICAL UTILITIES
// -----------------------------------------------------------------------------

// count digits in a number
int count_digits(int n) {
    return to_string(n).size();
};

// recursive power function.
double power(double num, int pow) {
    if (pow == 0) {
        return 1;
    }
    return num * power(num, pow - 1);
}

// factorial calculation.
int factorial(int num) {
    int sum = 1;
    for (int i = 0; i < num; ++i) {
        sum *= i+1;
    }
    return sum;
};

// combination calculation for given numbers.
double combination(int total_number, int pull_number) {
    return factorial(total_number) / (factorial(pull_number) * factorial(total_number - pull_number));
};

// computes binomial probability.
double binomial_probability(double probability, int total_number, int pull_number) {
    double result = 1.0;
    for (int i = 0; i < pull_number; ++i) {
        result *= (total_number - i);
        result /= (i + 1);
    }
    return result * pow(probability, pull_number) * pow(1 - probability, total_number - pull_number);
}

// computes the mean for a probability distribution
double mean(vector<double> &values, vector<double> &probabilities) {
    if (values.size() != probabilities.size()) {
        throw std::invalid_argument("not equal");
    }
    double result = 0;
    for (int i = 0; i < values.size(); ++i) {
        result += values[i] * probabilities[i];
    }
    return result;
}

// computes the variance for a porbability distribution
double variance(vector<double> values, vector<double> probabilities) {
    if (values.size() != probabilities.size()) {
        throw std::invalid_argument("not equal");
    }
    double result = 0;
    double m = mean(values, probabilities);
    for (int i = 0; i < values.size(); ++i) {
        result += power(values[i] - m, 2) * probabilities[i];
    }
    return result;
};

// computes the standard deviation of a probability distribution
double standard_deviation(vector<double> values, vector<double> probabilities) {
    return sqrt(variance(values, probabilities));
};

// -----------------------------------------------------------------------------
// 3. TEMPLATE UTILITIES
// -----------------------------------------------------------------------------

// count items in a vector.
template <typename T>
int how_many_items(const std::vector<T>& things) {
    int i = 0;
    for (T thing : things) {
        i++;
    }
    return i;
};
template <typename T, size_t N>
std::vector<T> array_to_vector(const T (&array)[N]) {
    std::vector<T> vector;
    for (int i = 0; i < N; i++) {
        vector.push_back(array[i]);
    }
    return vector;
};

// -----------------------------------------------------------------------------
// 4. TEXT HANDLING
// -----------------------------------------------------------------------------

// write data to a text file.
template <typename T>
void write_to_text(const T& data, const std::string& filename, std::ios_base::openmode mode = std::ios::app) {
    std::ofstream outfile;
    outfile.open(filename, mode);
    if (!outfile.is_open()) {
        throw std::runtime_error("failed to open the file");
    }
    outfile << data;
    outfile.close();
}

// read data from a text file.
std::string read_from_text(const std::string& filename, std::ios::openmode mode = std::ios::app) {
    std::ifstream infile;
    infile.open(filename, mode);
    if (!infile.is_open()) {
        throw std::runtime_error("failed to open the file");
    }
    std::string content((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>());
    infile.close();

    return content;
}

// XOR encryption for strings.
std::string xor_text_encryption(std::string& data, const std::string& key) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] ^= key[i % key.size()];
    }
    return data;
}

// -----------------------------------------------------------------------------
// 5. BINARY HANDLING
// -----------------------------------------------------------------------------

// write data to a binary file.
template <typename T>
void write_to_binary(T& data, const std::string& filename, std::ios::openmode mode = std::ios::trunc) {
    std::ofstream outfile;
    outfile.open(filename, mode);
    if (!outfile.is_open()) {
        throw std::runtime_error("failed to open the file");
    }
    outfile.write(reinterpret_cast<char*>(&data), sizeof(data));
    outfile.close();
}

// read data from a binary file.
std::vector<std::uint8_t> read_from_binary(const std::string& filename, std::ios::openmode mode = std::ios::trunc) {
    std::ifstream infile;
    infile.open(filename, std::ios::binary | mode);
    if (!infile.is_open()) {
        throw std::runtime_error("failed to open the file");
    }

    infile.seekg(0, std::ios::end);
    std::streamsize size = infile.tellg();
    infile.seekg(0, std::ios::beg);

    std::vector<uint8_t> buffer(size);

    if (!infile.read(reinterpret_cast<char*>(buffer.data()), size)) {
        throw std::runtime_error("failed to read the file");
    }

    infile.close();

    return buffer;
};

// XOR encryption for binary data.
std::vector<std::uint8_t> xor_binary_encryption(std::vector<std::uint8_t> data, const std::string& key) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] ^= key[i % key.size()];
    }
    return data;
}

#endif //USEFUL_USEFUL_H
