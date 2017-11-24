//
// Created by EL on 23.05.17.
//

#include <string>
#include <vector>
#include <boost/math/special_functions/gamma.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../PValue.h"
#include "../InputFile.h"


// Calculating the number of unique elements in a phenotype vector
TEST(pvalue_check, calculate_number_of_elements_vector_checked){
    vector<unsigned short> test_vec = {0, 2, 1, 1, 0, 2, 1, 0, 1, 0, 2, 1};
    PhenotypeStatistics phen_stat = PValue::calcNumElem(test_vec);
    ASSERT_EQ(3, phen_stat.num_elements);
};

// Calculating the number of unique elements in a genotype vector
TEST(pvalue_check, calculate_number_of_elements_genotype_vector_checked){
    vector<unsigned short> test_matrix = {1, 0, 1, 0, 2, 1};
    vector<unsigned short> test_vector = {3, 0, 0, 1, 1, 2};
    vector<vector<int>> V(4, vector<int>(4));
    PValue::fillVMatrix(test_matrix, test_vector, V);
    int result = PValue::calcNumElementsInGenotype(V);
    ASSERT_EQ(3, result);
};

// Calculating the number of unique elements in a phenotype vector filled with zeros
TEST(pvalue_check, calculate_number_of_elements_vector_all_zeros_checked){
    vector<unsigned short> test_vec = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    PhenotypeStatistics phen_stat = PValue::calcNumElem(test_vec);
    ASSERT_EQ(1, phen_stat.num_elements);
};

// Calculating the number of unique elements in a genotype vector filled with zeros
TEST(pvalue_check, calculate_number_of_elements_genotype_vector_all_zeros_checked){
    vector<unsigned short> test_matrix = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    vector<unsigned short> test_vector = {3, 0, 0, 1, 1, 2, 2, 2, 1, 0, 3, 0};
    vector<vector<int>> V(4, vector<int>(4));
    PValue::fillVMatrix(test_matrix, test_vector, V);
    int result = PValue::calcNumElementsInGenotype(V);
    ASSERT_EQ(1, result);
};

// Calculating the number of unique elements in a genotype vector where all elements are '3'
// The element '3' should not be considered
TEST(pvalue_check, calculate_number_of_elements_genotype_vector_all_threes_checked){
    vector<unsigned short> test_matrix = {3, 3, 3, 3, 3, 3};
    vector<unsigned short> test_vector = {3, 0, 0, 1, 1, 2};
    vector<vector<int>> V(4, vector<int>(4));
    PValue::fillVMatrix(test_matrix, test_vector, V);
    int result = PValue::calcNumElementsInGenotype(V);
    ASSERT_EQ(0, result);
};

// Calculating the number of unique elements in an empty vector
TEST(pvalue_check, calculate_number_of_elements_empty_vector_check){
    vector<unsigned short> empty = {};
    PhenotypeStatistics phen_stat = PValue::calcNumElem(empty);
    ASSERT_EQ(0, phen_stat.num_elements);
};

// Calculating the number of unique elements in an empty genotype vector
TEST(pvalue_check, calculate_number_of_elements_empty_genotype_vector_check){
    vector<unsigned short> empty_matrix = {};
    vector<unsigned short> empty_vector = {};
    vector<vector<int>> V(4, vector<int>(4));
    PValue::fillVMatrix(empty_matrix, empty_vector, V);
    int result = PValue::calcNumElementsInGenotype(V);
    ASSERT_EQ(0, result);
};

// Getting the value of the max element in phenotype vector
TEST(pvalue_check, calculate_max_element_in_phenotype_vector_check){
    vector<unsigned short> phenotype = {0,0,4,3,3,2,3,5,2,0,1,1,0};
    PhenotypeStatistics phen_stat = PValue::calcNumElem(phenotype);
    ASSERT_EQ(5, phen_stat.max_element);
};

// Getting the value of the max element in phenotype vector (max is last)
TEST(pvalue_check, calculate_max_element_last_in_phenotype_vector_check){
    vector<unsigned short> phenotype = {0,0,4,3,3,2,3,5,2,0,1,1,6};
    PhenotypeStatistics phen_stat = PValue::calcNumElem(phenotype);
    ASSERT_EQ(6, phen_stat.max_element);
};

// Getting the value of the max element in an empty phenotype vector
TEST(pvalue_check, calculate_max_element_in_empty_phenotype_vector_check){
    vector<unsigned short> phenotype = {};
    PhenotypeStatistics phen_stat = PValue::calcNumElem(phenotype);
    ASSERT_EQ(0, phen_stat.max_element);
};

// Checking the lower gamma function from the boost library
// The returned value for the given parameters is compared with the value obtained on this resource:
// http://keisan.casio.com/exec/system/1180573447
TEST(pvalue_check, lower_gamma_function_check){
    int p1 = 6;
    double p2 = 8.0;
    ASSERT_NEAR(97.0516, boost::math::tgamma_lower(p1, p2), 0.001);
};

// Testing the function prepareData for test type = "cd"
TEST(pvalue_check, prepare_data_cd_check){
    vector<vector<unsigned short>> cur_G = {};
    vector<unsigned short> cur_A = {};
    InputFile genotype("test.txt");
    vector<unsigned short> A = {1, 3, 0, 2, 1, 0};

    vector<vector<unsigned short>> result_G = {{0, 0, 0, 1, 1 ,0}, {2, 1, 0, 0, 2, 1},
                                              {1, 0, 0, 2, 2, 1}, {0, 0, 1, 1, 1, 0}};
    vector<unsigned short> result_A = {1, 3, 0, 2, 1, 0};

    TableEntry cur_test = {"","","cd",2,5,0};
    PValue::prepareData(cur_G, cur_A, genotype, A, cur_test);

    ASSERT_EQ(result_A, cur_A);
    ASSERT_EQ(result_G, cur_G);
}

// Testing the function prepareData for test type = "r"
TEST(pvalue_check, prepare_data_r_check){
    vector<vector<unsigned short>> cur_G = {};
    vector<unsigned short> cur_A = {};
    InputFile genotype("test.txt");
    vector<unsigned short> A = {1, 3, 0, 2, 1, 0};

    vector<vector<unsigned short>> result_G = {{0, 0, 0, 0, 0 ,0}, {1, 0, 0, 0, 1, 0},
                                              {0, 0, 0, 1, 1, 0}, {0, 0, 0, 0, 0, 0}};
    vector<unsigned short> result_A = {1, 3, 0, 2, 1, 0};

    TableEntry cur_test = {"","","r",2,5,0};
    PValue::prepareData(cur_G, cur_A, genotype, A, cur_test);

    ASSERT_EQ(result_A, cur_A);
    ASSERT_EQ(result_G, cur_G);
}

// Testing the function prepareData for test type = "d"
TEST(pvalue_check, prepare_data_d_check){
    vector<vector<unsigned short>> cur_G = {};
    vector<unsigned short> cur_A = {};
    InputFile genotype("test.txt");
    vector<unsigned short> A = {1, 3, 0, 2, 1, 0};

    vector<vector<unsigned short>> result_G = {{0, 0, 0, 1, 1 ,0}, {1, 1, 0, 0, 1, 1},
                                              {1, 0, 0, 1, 1, 1}, {0, 0, 1, 1, 1, 0}};
    vector<unsigned short> result_A = {1, 3, 0, 2, 1, 0};

    TableEntry cur_test = {"","","d",2,5,0};
    PValue::prepareData(cur_G, cur_A, genotype, A, cur_test);

    ASSERT_EQ(result_A, cur_A);
    ASSERT_EQ(result_G, cur_G);
}

// Testing the function prepareData for test type = "a"
TEST(pvalue_check, prepare_data_a_check){
    vector<vector<unsigned short>> cur_G = {};
    vector<unsigned short> cur_A = {};
    InputFile genotype("test.txt");
    vector<unsigned short> A = {1, 3, 0, 2, 1, 0};

    vector<vector<unsigned short>> result_G = {{0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1 ,0},
                                              {1, 0, 0, 0, 1, 0,   1, 1, 0, 0, 1, 1},
                                              {0, 0, 0, 1, 1, 0,   1, 0, 0, 1, 1, 1},
                                              {0, 0, 0, 0, 0, 0,   0, 0, 1, 1, 1, 0}};

    vector<unsigned short> result_A = {1, 3, 0, 2, 1, 0,  1, 3, 0, 2, 1, 0};

    TableEntry cur_test = {"","","a",2,5,0};
    PValue::prepareData(cur_G, cur_A, genotype, A, cur_test);

    ASSERT_EQ(result_A, cur_A);
    ASSERT_EQ(result_G, cur_G);
}

// Testing how a vector of strings is converted to a vector of unsigned chars
TEST(pvalue_check, map_phenotype_values_to_char_check){
    vector<string> phen = {"blue", "red", "green", "red", "black", "blue", "green"};
    vector<unsigned short> new_phen;
    new_phen = PValue::mapPhenotypeValuesToChar(phen);
    vector<unsigned short> result = {0, 1, 2, 1, 3, 0, 2};
    ASSERT_EQ(result, new_phen);
}

// Testing how an empty vector of strings is converted to a vector of unsigned chars
TEST(pvalue_check, map_phenotype_values_to_char_empty_vector_check) {
    vector<string> phen = {};
    vector<unsigned short> new_phen;
    new_phen = PValue::mapPhenotypeValuesToChar(phen);
    vector<unsigned short> result = {};
    ASSERT_EQ(result, new_phen);
}

// Testing fillVMatrix with 3x4 matrix
TEST(pvalue_check, fill_v_matrix_3x4_check){
    vector<vector<int>> V1(4, vector<int>(4)), V2(4, vector<int>(4)), V3(4, vector<int>(4)),
            V4(4, vector<int>(4)), real_V1, real_V2, real_V3, real_V4;
    vector<vector<unsigned short>> genotype = {{0, 0, 0, 1, 1, 0}, {2, 1, 0, 0, 2, 1},
                                              {1, 0, 0, 2, 2, 1}, {0, 0, 1, 1, 1, 0}};
    vector<unsigned short> phenotype = {1, 3, 0, 2, 1, 0};

    PValue::fillVMatrix(genotype[0], phenotype, V1);
    PValue::fillVMatrix(genotype[1], phenotype, V2);
    PValue::fillVMatrix(genotype[2], phenotype, V3);
    PValue::fillVMatrix(genotype[3], phenotype, V4);

    real_V1 = {{2, 1, 0, 1}, {0, 1, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    real_V2 = {{1, 0, 1, 0}, {1, 0, 0, 1}, {0, 2, 0, 0}, {0, 0, 0, 0}};
    real_V3 = {{1, 0, 0, 1}, {1, 1, 0, 0}, {0, 1, 1, 0}, {0, 0, 0, 0}};
    real_V4 = {{1, 1, 0, 1}, {1, 1, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

    ASSERT_EQ(real_V1, V1);
    ASSERT_EQ(real_V2, V2);
    ASSERT_EQ(real_V3, V3);
    ASSERT_EQ(real_V4, V4);
}

// Testing fillVMatrix with 2x4 matrix
TEST(pvalue_check, fill_v_matrix_2x4_check){
    vector<vector<int>> V1(4, vector<int>(4)), V2(4, vector<int>(4)), V3(4, vector<int>(4)),
             real_V1, real_V2, real_V3;
    vector<vector<unsigned short>> genotype = {{1, 1, 0, 0, 0, 1}, {0, 0, 0, 1, 1, 0},
                                              {1, 0, 0, 0, 0, 1}};
    vector<unsigned short> phenotype = {1, 3, 0, 2, 1, 0};

    PValue::fillVMatrix(genotype[0], phenotype, V1);
    PValue::fillVMatrix(genotype[1], phenotype, V2);
    PValue::fillVMatrix(genotype[2], phenotype, V3);

    real_V1 = {{1, 1, 1, 0}, {1, 1, 0, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    real_V2 = {{2, 1, 0, 1}, {0, 1, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    real_V3 = {{1, 1, 1, 1}, {1, 1, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

    ASSERT_EQ(real_V1, V1);
    ASSERT_EQ(real_V2, V2);
    ASSERT_EQ(real_V3, V3);
}

// Testing fillVMatrix with '3' elements in genotype
TEST(pvalue_check, fill_v_matrix_with_3_in_genotype_check){
    vector<vector<int>> V1(4, vector<int>(4)), V2(4, vector<int>(4)), V3(4, vector<int>(4)),
            V4(4, vector<int>(4)), real_V1, real_V2, real_V3, real_V4;
    vector<vector<unsigned short>> genotype = {{0, 3, 0, 1, 1, 0}, {2, 1, 3, 3, 2, 1},
                                              {1, 0, 0, 2, 3, 1}, {0, 0, 3, 1, 1, 0}};
    vector<unsigned short> phenotype = {2, 3, 0, 1, 1, 0};

    PValue::fillVMatrix(genotype[0], phenotype, V1);
    PValue::fillVMatrix(genotype[1], phenotype, V2);
    PValue::fillVMatrix(genotype[2], phenotype, V3);
    PValue::fillVMatrix(genotype[3], phenotype, V4);

    real_V1 = {{2, 0, 1, 0}, {0, 2, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 1}};
    real_V2 = {{0, 0, 0, 0}, {1, 0, 0, 1}, {0, 1, 1, 0}, {1, 1, 0, 0}};
    real_V3 = {{1, 0, 0, 1}, {1, 0, 1, 0}, {0, 1, 0, 0}, {0, 1, 0, 0}};
    real_V4 = {{1, 0, 1, 1}, {0, 2, 0, 0}, {0, 0, 0, 0}, {1, 0, 0, 0}};

    ASSERT_EQ(real_V1, V1);
    ASSERT_EQ(real_V2, V2);
    ASSERT_EQ(real_V3, V3);
    ASSERT_EQ(real_V4, V4);
}

// Testing calculateChiSqr when V is 2x4
TEST(pvalue_check, calculate_chi_sqr_V2x4){
    int V_rows = 2, V_cols = 4;
    vector<vector<int>> V;
    V = {{3, 4, 1, 0}, {1, 0, 1, 2}};

    double chi_sqr, real_chi_sqr = 6.375;
    chi_sqr = PValue::calculateChiSqr(V, V_rows, V_cols);

    ASSERT_NEAR(real_chi_sqr, chi_sqr, 0.001);
}

// Testing calculateChiSqr when V is 3x4
TEST(pvalue_check, calculate_chi_sqr_V3x4){
    int V_rows = 3, V_cols = 4;
    vector<vector<int>> V;
    V = {{1, 0, 0, 1}, {1, 0, 1, 0}, {0, 1, 0, 0}};

    double chi_sqr, real_chi_sqr = 7.5;
    chi_sqr = PValue::calculateChiSqr(V, V_rows, V_cols);

    ASSERT_NEAR(real_chi_sqr, chi_sqr, 0.001);
}

// Testing calculateChiSqr when V is 3x4 and the rows are proportional
TEST(pvalue_check, calculate_chi_sqr_V3x4_proportional){
    int V_rows = 3, V_cols = 4;
    vector<vector<int>> V;
    V = {{1, 1, 0, 0}, {2, 2, 0, 0}, {1, 1, 0, 0}};

    double chi_sqr, real_chi_sqr = 0.0;
    chi_sqr = PValue::calculateChiSqr(V, V_rows, V_cols);

    ASSERT_NEAR(real_chi_sqr, chi_sqr, 0.001);
}

// Testing calculateChiSqr when V is empty
TEST(pvalue_check, calculate_chi_sqr_empty_matrix){
    int V_rows = 0, V_cols = 0;
    vector<vector<int>> V;
    V = {};

    double chi_sqr, real_chi_sqr = 0.0;
    chi_sqr = PValue::calculateChiSqr(V, V_rows, V_cols);

    ASSERT_NEAR(real_chi_sqr, chi_sqr, 0.001);
}

// Testing calculateChiSqr when V has a row of zeros
TEST(pvalue_check, calculate_chi_sqr_row_of_zeros){
    int V_rows = 4, V_cols = 4;
    vector<vector<int>> V;
    V = {{1, 0, 0, 1}, {1, 0, 1, 0}, {0, 0, 0, 0}, {0, 1, 0, 0}};

    double chi_sqr, real_chi_sqr = 7.5;
    chi_sqr = PValue::calculateChiSqr(V, V_rows, V_cols);

    ASSERT_NEAR(real_chi_sqr, chi_sqr, 0.001);
}

// Testing calculateChiSqr when V has a column of zeros
TEST(pvalue_check, calculate_chi_sqr_column_of_zeros){
    int V_rows = 3, V_cols = 5;
    vector<vector<int>> V;
    V = {{1, 0, 0, 0, 1}, {1, 0, 1, 0, 0}, {0, 1, 0, 0, 0}};

    double chi_sqr, real_chi_sqr = 7.5;
    chi_sqr = PValue::calculateChiSqr(V, V_rows, V_cols);

    ASSERT_NEAR(real_chi_sqr, chi_sqr, 0.001);
}

// Testing calcPValue
TEST(pvalue_check, calc_p_value_check){
    vector<vector<unsigned short>> genotype = {{0, 1, 0, 0, 0, 0,   0, 1, 1, 1, 0, 0},
                                              {0, 1, 0, 0, 1, 0,   0, 1, 1, 0, 1, 1},
                                              {0, 0, 0, 1, 0, 0,   0, 1, 1, 1, 1, 1}};

    vector<unsigned short> phenotype = {1, 3, 0, 2, 1, 0,  1, 3, 0, 2, 1, 0};
    double pValue = PValue::calcPValue(genotype, phenotype, "a");
    double realPValue = 2.2095956;

    ASSERT_NEAR(realPValue, pValue, 0.0000001);
}

// Testing calcPValue when there is no '0'-elements in the phenotype
TEST(pvalue_check, calc_p_value_cut_phenotype_check){
    vector<vector<unsigned short>> genotype = {{0, 1, 0, 0, 0, 0,   0, 1, 1, 1, 0, 0},
                                               {0, 1, 0, 0, 1, 0,   0, 1, 1, 0, 1, 1},
                                               {0, 0, 0, 1, 0, 0,   0, 1, 1, 1, 1, 1}};

    vector<unsigned short> phenotype = {2, 4, 1, 3, 2, 1,  2, 4, 1, 3, 2, 1};
    double pValue = PValue::calcPValue(genotype, phenotype, "a");
    double realPValue = 2.2095956;

    ASSERT_NEAR(realPValue, pValue, 0.0000001);
}

/*
//  Main test
TEST(pvalue_check, adjust_p_value_check){
    vector<TestsData> tests;
    tests = {{"cd",2,5}, {"r",6,8}, {"d",10,11}, {"a",12,14}, {"d",16,19}};
    ExecutionParameters parameters;
    parameters.maxReplications = 10; // Should be around 10^(-8)
    parameters.k = 10;
    parameters.isAdaptive = true;
    InputFile genotype("test.txt");
    vector<unsigned short> phenotype = {1, 3, 0, 2, 1, 0};

    vector<double> results = PValue::adjustPValue(tests, genotype, phenotype, parameters);
    vector<double> real_results = {0.8, 0.2, 0, 0.9, 0.9};

    for (int j = 0; j<real_results.size(); ++j) {
        ASSERT_NEAR(real_results[j], results[j], 0.001);
    }
}
 */

// Checking error message when the given file is empty.
TEST(pvalue_check, prepare_data_check_empty_genotype_matrix){
    InputFile genotype("emptyfile.txt");
    vector<vector<unsigned short>> cur_G;
    vector<unsigned short> cur_A;
    vector<unsigned short> const A = {1,1,3,2,4,5};
    TableEntry cur_test = {"","","cd", 0, 24,0};

    ASSERT_DEATH(PValue::prepareData(cur_G, cur_A, genotype, A, cur_test), "Genotype file is empty!");
}

// Checking error message when the given phenotype is empty.
TEST(pvalue_check, prepare_data_check_empty_phenotype_vector){
    InputFile genotype("test.txt");
    vector<vector<unsigned short>> cur_G;
    vector<unsigned short> cur_A;
    vector<unsigned short> const A = {};
    TableEntry cur_test = {"","","cd", 0, 24,0};

    ASSERT_DEATH(PValue::prepareData(cur_G, cur_A, genotype, A, cur_test), "Phenotype is empty!");
}

// Checking the gamma function from the standard library
// The returned value for the given parameters is compared with the value returned by gamma() function
// in MATLAB
TEST(pvalue_check, gamma_function_check){
    int p = 6;
    ASSERT_NEAR(120, tgamma(p), 0.001);
};

TEST(pvalue_check, check_num_elem){
    vector<unsigned short> const & genotype = {0, 1, 1, 2, 1, 3};
    bool answ = PValue::checkNumElem(genotype, 3);
    bool real_answ = true;
    ASSERT_EQ(real_answ, answ);
}

TEST(pvalue_check, check_num_elem_2){
    vector<unsigned short> const & genotype = {3, 1, 1, 2, 1, 2};
    bool answ = PValue::checkNumElem(genotype, 3);
    bool real_answ = true;
    ASSERT_EQ(real_answ, answ);
}

TEST(pvalue_check, check_num_elem_3){
    vector<unsigned short> const & genotype = {1, 1, 1, 3, 1, 3};
    bool answ = PValue::checkNumElem(genotype, 3);
    bool real_answ = false;
    ASSERT_EQ(real_answ, answ);
}

TEST(pvalue_check, check_num_elem_4){
    vector<unsigned short> const & genotype = {1, 1, 1, 3, 1, 1};
    bool answ = PValue::checkNumElem(genotype, 3);
    bool real_answ = false;
    ASSERT_EQ(real_answ, answ);
}

TEST(pvalue_check, check_num_elem_5){
    vector<unsigned short> const & genotype = {1, 1, 1, 1, 1, 1};
    bool answ = PValue::checkNumElem(genotype, 3);
    bool real_answ = false;
    ASSERT_EQ(real_answ, answ);
}