//
// Created by doby on 17/11/17.
//

#ifndef APV_C_VERSION_CALCPVALUE_H
#define APV_C_VERSION_CALCPVALUE_H

//#include "InputData.h"
//#include <gtest/gtest.h>
//#include <gmock/gmock.h>
//#include <iostream>
#include "stdio.h"
#include <stdbool.h>
#include <math.h>


// Defines the grouping of the genotype = {"cd", "d", "r", "a"}
typedef enum AlternativeHypothesisType {
    eCD, // Codominant
    eR,  // Dominant
    eD,  // Recessive
    eA   // Allelic
} AlternativeHypothesisType;

// New format for the Top Hits Table
typedef struct TableEntry {
    char *test_index;
    char *chromosome_index;
    char *ID;                 // Defines the grouping of the genotype. ID = {"cd", "d", "r", "a"}
    size_t lower;
    size_t upper;
    double adjusted_p_value;  // Place for the results of computations
} TableEntry;

typedef struct PhenotypeStatistics {
    unsigned short max_element;
    int num_elements;
} PhenotypeStatistics;

// Execution parameters
typedef struct ExecutionParameters {
    int maxReplications; // Maximum amount of replications
    int k; // Maximum amount of P-values (D_cur) that are larger than the initial P-value (D_main)
    bool isAdaptive; // Defines whether K will be used
} ExecutionParameters;

//void adjustPValue(int type);
void adjustPValue(TableEntry *tests, const unsigned short *A, size_t tests_len, size_t A_len, const ExecutionParameters cont, char *path, int type); // used to contain InputData &G, can be substituted by a pointer to a function
void prepareData(unsigned short **cur_G, unsigned short **cur_A, const unsigned short *A, size_t *G_len, size_t *A_len, TableEntry cur_test, char *path, int type); // InputData & G
double calcPValue(const unsigned short *cur_G, const unsigned short *cur_A, char *ID, size_t G_len, size_t A_len);
unsigned short *createGenotypeMatrix(char *path, size_t G_len, size_t lower, size_t upper, int type);
bool checkNumElemInGenotype(const unsigned short *genotype, size_t G_len);
AlternativeHypothesisType hashIt (const char *inString);
PhenotypeStatistics calcNumElem(const unsigned short *phenotype, size_t phen_len);
unsigned short* doubleSizeOfPhenotype(unsigned short *cur_A, size_t A_len);
unsigned short* doubleSizeOfGenotype(unsigned short *cur_G, size_t G_len, size_t row_num);
PhenotypeStatistics calcNumElem(const unsigned short *phenotype, size_t phen_len);
void fillVMatrix(const unsigned short *cur_G, const unsigned short *cur_A, int *V, int V_rows, int V_cols, size_t col_num);
double calculateChiSqr(const int *V, int V_rows, int V_cols);
void phenotypeRandomPermutation(unsigned short *phenotype, size_t phenotypeLength);
int calcNumElementsInGenotype(const int *V, int rowNum, int colNum); // maybe 2d array
size_t rand_interval(size_t min, size_t max);



#endif //APV_C_VERSION_CALCPVALUE_H
