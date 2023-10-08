#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include "baseline_dtw.hpp"
#include "dtw.hpp"
#include <assert.h>

using namespace std;

vector<double> generate_random_vector(int num_elements, unsigned int seed) {
    // Create a Mersenne Twister random number generator with the given seed
    std::mt19937 rng(seed);

    // Create a uniform distribution for generating random doubles in the range [-2.5, 2.5)
    std::uniform_real_distribution<double> dist(-2.5, 2.5);

    // Generate the random values and store them in a vector
    std::vector<double> random_vector(num_elements);
    for (int i = 0; i < num_elements; i++) {
        random_vector[i] = dist(rng);
    }

    return random_vector;
}

double baseline_dtw(vector<double> a, vector<double> b){
    vector<vector<double>> internal_a;
    for(int i = 0; i < a.size(); i++){
        internal_a.push_back({a[i]});
    }

    vector<vector<double>> internal_b;
    for(int i = 0; i < b.size(); i++){
        internal_b.push_back({b[i]});
    }

    return DTW::dtw_distance_only(internal_a, internal_b, 1);
}

std::vector<float> convert_to_float_vector(const std::vector<double>& double_vector) {
    std::vector<float> float_vector(double_vector.size());

    // Use std::transform to apply a lambda function that converts each double element to a float
    std::transform(double_vector.begin(), double_vector.end(), float_vector.begin(), [](double d) {
        return static_cast<float>(d);
    });

    return float_vector;
}

double rawalign_dtw(vector<double> a, vector<double> b,
    float f(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element)){
    vector<float> float_a = convert_to_float_vector(a);
    vector<float> float_b = convert_to_float_vector(b);
    return f(float_a.data(), float_a.size(), float_b.data(), float_b.size(), false);
}

double rawalign_dtw_banded(vector<double> a, vector<double> b, uint32_t band_radius,
    float f(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, int band_radius, bool exclude_last_element)){
    vector<float> float_a = convert_to_float_vector(a);
    vector<float> float_b = convert_to_float_vector(b);
    return f(float_a.data(), float_a.size(), float_b.data(), float_b.size(), band_radius, false);
}

dtw_result rawalign_dtw_tb(vector<double> a, vector<double> b,
    dtw_result f(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element)){
    vector<float> float_a = convert_to_float_vector(a);
    vector<float> float_b = convert_to_float_vector(b);
    return f(&float_a[0], float_a.size(), &float_b[0], float_b.size(), false);
}

void run_comparison(){
    unsigned int seed = 42;
    vector<double> a = generate_random_vector(10, seed);
    vector<double> b = generate_random_vector(10, seed+1);
    //vector<double> a = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    //vector<double> b = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    
    cout << "baseline_dtw                              " << baseline_dtw(a, b) << endl;
    cout << "DTW_global                                " << rawalign_dtw(a, b, DTW_global) << endl;
    cout << "DTW_global_slow                           " << rawalign_dtw(a, b, DTW_global_slow) << endl;
    cout << "DTW_global_diagonalbanded                 " << rawalign_dtw_banded(a, b, 3, DTW_global_diagonalbanded) << endl;
    cout << "DTW_global_slantedbanded                  " << rawalign_dtw_banded(a, b, 3, DTW_global_slantedbanded) << endl;
    cout << "DTW_global_slantedbanded_antidiagonalwise " << rawalign_dtw_banded(a, b, 3, DTW_global_slantedbanded_antidiagonalwise) << endl;

    cout << endl;

    cout << "DTW_semiglobal      " << rawalign_dtw(a, b, DTW_semiglobal) << endl;
    cout << "DTW_semiglobal_slow " << rawalign_dtw(a, b, DTW_semiglobal_slow) << endl;

    cout << endl;

    {
        dtw_result result = rawalign_dtw_tb(a, b, DTW_global_tb);
        cout << "DTW_global_tb   " << result.cost << endl;
    }
    {
        dtw_result result = rawalign_dtw_tb(a, b, DTW_semiglobal_tb);
        cout << "DTW_semiglobal_tb   " << result.cost << endl;
    }
}

void run_comparison_tb(){
    unsigned int seed = 42;
    vector<double> a = generate_random_vector(10, seed);
    vector<double> b = generate_random_vector(10, seed+1);

    {
        dtw_result result = rawalign_dtw_tb(a, b, DTW_global_tb);
        cout << "DTW_global_tb   " << result.cost << endl;
        for(auto& element : result.alignment){
            cout << element.position.i << " " << element.position.j << " " << element.difference << endl;
        }
    }

    {
        dtw_result result = rawalign_dtw_tb(a, b, DTW_semiglobal_tb);
        cout << "DTW_global_tb   " << result.cost << endl;
        for(auto& element : result.alignment){
            cout << element.position.i << " " << element.position.j << " " << element.difference << endl;
        }
    }
}

int get_necessary_band_radius(dtw_result aln){
    float target_slope = aln.alignment.back().position.j / (float)aln.alignment.back().position.i;
    int max_diff = 0;
    for(int i = 0; i < aln.alignment.size(); i++){
        int diff_from_center_slope = (int)ceil(abs(aln.alignment[i].position.j - aln.alignment[i].position.i * target_slope));
        max_diff = max(max_diff, diff_from_center_slope);
    }
    return max_diff;
}

#define APPROX_EQ(a, b) (abs(a - b) < 0.001)
bool random_unit_test(uint32_t a_length, uint32_t b_length, uint32_t seed){
    vector<double> a = generate_random_vector(a_length, seed);
    vector<double> b = generate_random_vector(b_length, seed+1);

    double baseline = baseline_dtw(a, b);
    bool passed = true;

    dtw_result dtw_tb = rawalign_dtw_tb(a, b, DTW_global_tb);
    if(!APPROX_EQ(dtw_tb.cost, baseline)){	
        passed = false;
    }
    int band_radius = max(0, get_necessary_band_radius(dtw_tb));

    vector<double> dtw_results;
    dtw_results.push_back(rawalign_dtw(a, b, DTW_global));
    dtw_results.push_back(rawalign_dtw(a, b, DTW_global_slow));
    //dtw_results.push_back(rawalign_dtw_banded(a, b, band_radius, DTW_global_diagonalbanded));
    dtw_results.push_back(rawalign_dtw_banded(a, b, band_radius, DTW_global_slantedbanded));
    dtw_results.push_back(rawalign_dtw_banded(a, b, band_radius, DTW_global_slantedbanded_antidiagonalwise));

    for(auto& result : dtw_results){
        if(!APPROX_EQ(result, baseline)){
            passed = false;
        }
    }

    if(!passed){
        //set high precision for cout
        cout << fixed << setprecision(6);
        cout << "***** FAILED *****" << endl;
        cout << "a_length: " << a_length << endl;
        cout << "b_length: " << b_length << endl;
        cout << "seed: " << seed << endl;
        cout << "baseline: " << baseline << endl;
        cout << "dtw_tb: " << dtw_tb.cost << endl;
        cout << "band_radius: " << band_radius << endl;
        for(auto& result : dtw_results){
            cout << "result: " << result << endl;
        }
    }

    return passed;
}

bool run_various_random_tests(int n_tests){
    int test_groups = 7;
    int tests_per_group = n_tests / test_groups;

    //tiny square tests
    for(int i = 0; i < tests_per_group; i++){
        if(!random_unit_test(4, 4, i)){
            return false;
        }
    }

    //small square tests
    for(int i = 0; i < tests_per_group; i++){
        if(!random_unit_test(10, 10, i)){
            return false;
        }
    }

    //small rect tests
    for(int i = 0; i < tests_per_group; i++){
        if(!random_unit_test(20, 10, i)){
            return false;
        }
    }

    //small irregular rect tests
    for(int i = 0; i < tests_per_group; i++){
        if(!random_unit_test(25, 10, i)){
            return false;
        }
    }

    //large square test
    for(int i = 0; i < tests_per_group; i++){
        if(!random_unit_test(100, 100, i)){
            return false;
        }
    }

    //large rect test
    for(int i = 0; i < tests_per_group; i++){
        if(!random_unit_test(200, 50, i)){
            return false;
        }
    }

    //large irregular rect test
    for(int i = 0; i < tests_per_group; i++){
        if(!random_unit_test(200, 30, i)){
            return false;
        }
    }

    return true;
}

double opt_blocker;
#define BENCHMARK(REPETITIONS, fcall) {\
    auto start = chrono::high_resolution_clock::now(); \
    opt_blocker = 0; \
    for(int i = 0; i < REPETITIONS; i++){ \
        opt_blocker += fcall; \
    } \
    auto end = chrono::high_resolution_clock::now(); \
    cout << #fcall << " time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() / (double)REPETITIONS << "us" << endl; \
}

void performance_benchmark(int iterations, int a_length, int b_length, float band_radius_fraction){
    unsigned int seed = 42;
    vector<double> a = generate_random_vector(a_length, seed);
    vector<double> b = generate_random_vector(b_length, seed+1);
    int band_radius = (int)(a.size() * band_radius_fraction);

    BENCHMARK(iterations, rawalign_dtw(a, b, DTW_global));
    BENCHMARK(iterations, rawalign_dtw(a, b, DTW_global_slow));
    BENCHMARK(iterations, rawalign_dtw_banded(a, b, band_radius, DTW_global_diagonalbanded));
    BENCHMARK(iterations, rawalign_dtw_banded(a, b, band_radius, DTW_global_slantedbanded));
    BENCHMARK(iterations, rawalign_dtw_banded(a, b, band_radius, DTW_global_slantedbanded_antidiagonalwise));
}

int main(int argc, char* argv[]){
    if(argc >= 2 && string(argv[1]) == "--performance-benchmark"){
        assert(argc == 6);
        int iterations = atoi(argv[2]);
        int a_length = atoi(argv[3]);
        int b_length = atoi(argv[4]);
        float band_radius_fraction = atof(argv[5]);
        performance_benchmark(iterations, a_length, b_length, band_radius_fraction);
        return 0;
    }

    //run_comparison();

    int n_tests = 10000;
    bool passed = run_various_random_tests(n_tests);
    if(!passed){
        return 1;
    }
    std::cout << "Passed " << n_tests << " random unit tests" << std::endl;

    //random_unit_test(25, 10, 24);

    //unsigned int seed = 42;
    //vector<double> a = generate_random_vector(10, seed);
    //vector<double> b = generate_random_vector(10, seed+1);
    //vector<double> a = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
    //vector<double> b = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
    //rawalign_dtw(a, b, DTW_global_SIMD16);
    //std::cout << rawalign_dtw(a, b, DTW_global_SIMD16) << std::endl;
    //std::cout << rawalign_dtw_banded(a, b, 2, DTW_global_slantedbanded) << std::endl;

    return 0;
}
