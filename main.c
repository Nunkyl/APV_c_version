#include <stdio.h>
#include <stdlib.h>
#include "PValue.h"
#include <gsl/gsl_sf_gamma.h>

#include <string.h>
#include <time.h>
#include <stddef.h>

bool checkNumElemInGenotype2_n(const uint8_t *genotype, size_t phn_len) {
    unsigned short p = 3;
    size_t elements_size = 4;
    unsigned short elements[4] = {0};
    unsigned short bits[4];
    int count = 0;
    union gen cur_gen;
    int bits_num = 4;

    for (size_t j = 0; j < (phn_len + 3)/4; ++j) {
        cur_gen.bin = genotype[j];
        bits[0] = cur_gen.n0, bits[1] = cur_gen.n1, bits[2] = cur_gen.n2, bits[3] = cur_gen.n3;
        if(j == (phn_len + 3)/4 - 1  && phn_len%4 != 0) bits_num = phn_len%4 + 1;
        for (int i = 0; i < bits_num; i++) {
            if (bits[i] != p) {
                if (elements[bits[i]] == 0) {
                    elements[bits[i]]++;
                    count++;
                    if (count == 2) {
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

PhenotypeStatistics calcNumElem_n(const unsigned short *phenotype, size_t phen_len){
    PhenotypeStatistics phen_stat = {.max_element = 0, .num_elements = 0};
    int elements_size = 10;
    unsigned short *elements = malloc(elements_size * sizeof(*elements));
    if (!elements) { /* Error handling */ }
    int count = 0;
    bool flag = 0;
    for(size_t j = 0; j < phen_len; j++) {
        flag = 0;
        for (int i = 0; i < count; ++i) {
            if (elements[i] == phenotype[j]) {flag = 1; break;}
        }
        if (!flag) {
            if (count == elements_size){
                elements_size = elements_size * 2;
                realloc(elements, elements_size * sizeof(*elements));
                if (!elements) { /* Error handling */  }
            }
            elements[count] = phenotype[j];
            count++;
            if (phenotype[j] > phen_stat.max_element) phen_stat.max_element = phenotype[j];
        }
    }
    free(elements);
    phen_stat.num_elements = count;
    return phen_stat;
}

void fillVMatrix_n(const uint8_t *cur_G, const unsigned short *cur_A, int *V, int V_rows, int V_cols, size_t col_num){  // col_num - The amount of patients

    union gen cur_gen;
    uint8_t bits[4];
    int bits_num = 4;

    // Fill V with zeros
    memset (V,0,V_rows * V_cols * sizeof(*V));

    for (int k = 0; k < V_rows * V_cols; ++k) {
        V[k] = 0;
    }

    // Fill V (G_car x A_car)
    for (size_t j = 0; j < (col_num + 3)/4; ++j) {
        cur_gen.bin = cur_G[j]; // uint8_t a = cur_G[j];
        bits[0] = cur_gen.n0, bits[1] = cur_gen.n1, bits[2] = cur_gen.n2, bits[3] = cur_gen.n3;
        if(j == (col_num + 3)/4 - 1 && col_num%4 != 0) bits_num = col_num%4 + 1;
        for (int i = 0; i < bits_num; i++) {
            V[bits[i] * V_cols + cur_A[j*4 + i]]++; // bits => a & 3; a >>= 2
        }
    }
}

double calculateChiSqr_n(const int *V, int V_rows, int V_cols){

    double chi_sqr = 0.0;
    int elem_num; // n in X^2
    int row_sum, col_sum;

    // Calc elem_num
    elem_num = 0;
    for (int m = 0; m < V_rows * V_cols; ++m) {
        elem_num += V[m];
    }
    for (int m = 0; m < V_rows; ++m) {
        row_sum = 0;
        for (int j = 0; j < V_cols; ++j) {
            row_sum = row_sum + V[m*V_cols + j];
        }
        for (int n = 0; n < V_cols; ++n) {
            col_sum = 0;
            for (int j = 0; j < V_rows; ++j) {
                col_sum = col_sum + V[j*V_cols + n];
            }
            if (row_sum*col_sum != 0) {
                double tmp = (V[m*V_cols + n] - ((double)(row_sum*col_sum))/elem_num);
                chi_sqr = chi_sqr + tmp * tmp / ((double)(row_sum*col_sum)/elem_num);
            }
        }
    }
    return chi_sqr;
}

// Calculates the exact number of elements in a genotype vector based on the V matrix
int calcNumElementsInGenotype_n(const int *V, int rowNum, int colNum){
    int num_elem = 0;
    int sum = 0;
    for (int i = 0; i < rowNum - 1; i++) {
        sum = 0;
        for (int j = 0; j < colNum; j++) {
            sum += V[i*colNum + j];
            if (sum > 0) {num_elem++; break;}
        }
    }
    return num_elem;
}

double calcPValue_n(const uint8_t *cur_G, const unsigned short *cur_A, size_t G_len, size_t A_len){

    double D = 0.0;     // Return value
    int D_num = 0;      // Amount of valid D values
    int V_rows, V_cols; // Size of V matrix
    int G_car;          // The cardinality of the sets of values of the phenotype and genotype
    int s;              // Parameters for calculating the gamma function
    double chi_sqr = 0.0;
    size_t row_num = G_len/A_len; // The amount of rows in the current genotype matrix
    PhenotypeStatistics phen_stat;
    size_t row_len_bytes = (A_len + 3)/4;

    size_t snp_len = (A_len + 3)/4;// Length of line in genotype matrix

    phen_stat = calcNumElem_n(cur_A, A_len);
    if (phen_stat.num_elements <= 1) return NAN; //?????????? quite NaN

    //Allocate memory for V matrix here, it will always be 4x(phen_stat.max_element+1)
    int *V = malloc(4*(phen_stat.max_element+1)*sizeof(*V));
    if (!V) { /* Error handling */  }

    for (int i = 0; i < row_num; ++i) {

        if (!checkNumElemInGenotype2_n(cur_G + i*(A_len+3)/4, A_len)) continue;
        fillVMatrix_n(cur_G + i*row_len_bytes,  cur_A, V, 4, (phen_stat.max_element+1), A_len);
        G_car = calcNumElementsInGenotype_n(V, 4, (phen_stat.max_element+1)); // Check this!

        // Calculate s and X^2
        s = (phen_stat.num_elements - 1)*(G_car - 1);
        chi_sqr = calculateChiSqr(V, 3, phen_stat.max_element+1); // Last line corresponding to '3' in genotype not considered

        D += -log10(gsl_sf_gamma_inc_Q(s, 2*chi_sqr)); // gsl_sf_gamma_inc(a,x)
        D_num++;
    }
    if (D_num == 0 || isnan(D)) return NAN; //?????????? quite NaN
    D = D/D_num;

    free(V);
    return D;
}

unsigned short * mapPhenotypeValuesToShort_n(const char *phn_name, size_t phn_cnt, ptrdiff_t phn_off){
    unsigned short *new_phn = malloc(phn_cnt *sizeof(*new_phn));
    if (!new_phn) { /* Error handling */  }
    size_t elements_cnt = 20, cur_num_elem = 0;
    char *elements = malloc(elements_cnt * phn_off * sizeof(*elements));
    if (!elements) { /* Error handling */  }
    bool flag = false;

    for(size_t i = 0; i < phn_cnt; i++){
        flag = false;
        for(size_t j = 0; j < cur_num_elem; j++){
            if(strncmp(phn_name + i*phn_off, elements + j*phn_off, phn_off) == 0){
                new_phn[i] = (unsigned short)j;
                flag = true;
            }
        }
        if(!flag){
            // Add memory if required
            if(cur_num_elem == elements_cnt){
                elements_cnt *= 2;
                realloc(elements, elements_cnt * sizeof(*elements));
                if (!elements) { /* Error handling */  }
            }

            strncpy(elements + cur_num_elem*phn_off, phn_name + i*phn_off, phn_off);
            new_phn[i] = (unsigned short)cur_num_elem;
            cur_num_elem++;
        }
    }
    free(elements);
    return new_phn;
}

void createGenotypeMatrix_uint8(uint8_t *genotype, char *path, size_t G_len, size_t lower, size_t upper, int type){
    if (type == 0) {  // If data comes from DB
        //gen = malloc (G_len * sizeof(*gen));
        //if (!gen) { /* Error handling */ return NULL; }

        union gen gen;
        gen.n0 = 0; gen.n1 = 0; gen.n2 = 0; gen.n3 = 1;
        genotype[0] = gen.bin;
        gen.n0 = 1; gen.n1 = 0; gen.n2 = 0; gen.n3 = 0;
        genotype[1] = gen.bin;

        gen.n0 = 2; gen.n1 = 1; gen.n2 = 0; gen.n3 = 0;
        genotype[2] = gen.bin;
        gen.n0 = 2; gen.n1 = 1; gen.n2 = 0; gen.n3 = 0;
        genotype[3] = gen.bin;

        gen.n0 = 1; gen.n1 = 0; gen.n2 = 0; gen.n3 = 2;
        genotype[4] = gen.bin;
        gen.n0 = 2; gen.n1 = 1; gen.n2 = 0; gen.n3 = 0;
        genotype[5] = gen.bin;

        gen.n0 = 0; gen.n1 = 0; gen.n2 = 1; gen.n3 = 1;
        genotype[6] = gen.bin;
        gen.n0 = 1; gen.n1 = 0; gen.n2 = 0; gen.n3 = 0;
        genotype[7] = gen.bin;
    }
    if (type == 1) {} // If data comes from text file
    if (type == 2) {} // If data comes from PLINK
    //return gen;
}

void prepareData_uint8(uint8_t **cur_G, unsigned short *cur_A, const unsigned short *A, size_t G_len, size_t A_len, TableEntry cur_test, char *path, int type){

    union gen cur_gen;
    unsigned short bits[4];

    size_t row_num = cur_test.upper - cur_test.lower + 1;
    createGenotypeMatrix_uint8(*cur_G, path, G_len, cur_test.lower, cur_test.upper, type); // Read the appropriate part of the genotype

    memcpy(cur_A, A, A_len * sizeof(cur_A));

    if (cur_test.ID == "a"){
        memcpy(cur_A + A_len/2, cur_A, A_len/2);
        uint8_t *buf_G =  malloc(G_len * sizeof(*buf_G)); // G_len = row_num * A_len
        if (!buf_G) { }

        size_t row_len = (A_len + 3)/4;
        for (size_t i = 0; i < row_num; i++){
            memcpy(buf_G + row_len*i*2, cur_G + row_len*i, row_len);
            memcpy(buf_G + row_len*(i*2+1), cur_G + row_len*i, row_len);
        }
        free(*cur_G);
        *cur_G = buf_G;
    }

    // Adjust the values of the genotype
    switch (hashIt(cur_test.ID)) {
        case eCD:
            break;
        case eR:
            for (size_t i = 0; i < G_len; i++) {
                cur_gen.bin = (*cur_G)[i];
                bits[0] = cur_gen.n0; bits[1] = cur_gen.n1; bits[2] = cur_gen.n2; bits[3] = cur_gen.n3;
                for(size_t j = 0; j < 4; j++) {
                    if (bits[j] == 1) bits[j] = 0;
                    if (bits[j] == 2) bits[j] = 1;
                }
                (*cur_G)[i] = cur_gen.bin;
            }
            break;
        case eD:
            for (size_t i = 0; i < G_len; i++) {
                cur_gen.bin = (*cur_G)[i];
                bits[0] = cur_gen.n0; bits[1] = cur_gen.n1; bits[2] = cur_gen.n2; bits[3] = cur_gen.n3;
                for(size_t j = 0; j < 4; j++) {
                    if (bits[j] == 2) bits[j] = 1;
                }
                (*cur_G)[i] = cur_gen.bin;
            }
            break;
        case eA:
            for (size_t i = 0; i < row_num; i++) {
                for (size_t j = 0; j < A_len/2; j++ ) {
                    if ((*cur_G)[i*A_len + j] == 1) (*cur_G)[i*A_len + j] = 0;
                    if ((*cur_G)[i*A_len + j] == 2) (*cur_G)[i*A_len + j] = 1;
                }
                for (size_t j = A_len/2; j < A_len; j++ ) {
                    if ((*cur_G)[i*A_len + j] == 2) (*cur_G)[i*A_len + j] = 1;
                }
            }

            size_t row_len = (A_len + 3)/4; // A_len = full length of one line in genotype matrix
            for (size_t i = 0; i < row_num; i++){
                for (size_t j = 0; j < row_len/2; j++ ) {
                    cur_gen.bin = (*cur_G)[i*row_len + j];
                    bits[0] = cur_gen.n0; bits[1] = cur_gen.n1; bits[2] = cur_gen.n2; bits[3] = cur_gen.n3;
                    for(size_t l = 0; l < 4; l++) {
                        if (bits[j] == 1) bits[j] = 0;
                        if (bits[j] == 2) bits[j] = 1;
                    }
                    (*cur_G)[i*row_len + j] = cur_gen.bin;
                }
                for (size_t j = row_len/2; j < row_len; j++ ){
                    cur_gen.bin = (*cur_G)[i*row_len + j];
                    bits[0] = cur_gen.n0; bits[1] = cur_gen.n1; bits[2] = cur_gen.n2; bits[3] = cur_gen.n3;
                    for(size_t l = 0; l < 4; l++) {
                        if (bits[j] == 2) bits[j] = 1;
                    }
                    (*cur_G)[i*row_len + j] = cur_gen.bin;
                }
            }

            break;
    }
}








int main() {


    // Test prepareData_uint8
    /*
    uint8_t *cur_G = malloc(8 * sizeof(*cur_G));
    unsigned short *cur_A = malloc(6 * sizeof(*cur_A));
    unsigned short *A = malloc(6 * sizeof(*A));
    A[0] = 1; A[1] = 3; A[2] = 0; A[3] = 2; A[4] = 1; A[5] = 0;
    size_t G_len = 8; // in bytes
    size_t A_len = 6; // in elements
    TableEntry cur_test; // Just one entry for now
    //tests[0] = {.test_index = "1", .chromosome_index = "1", .ID = "cd", .lower = 0, .upper = 3, .adjusted_p_value = 0.0};
    cur_test.test_index = "1";
    cur_test.chromosome_index = "1";
    cur_test.ID = "cd"; // r, d, a, d results: 0.7 1.0 0.6 0.2 0.0
    cur_test.lower = 0;
    cur_test.upper = 3;
    cur_test.adjusted_p_value = 0.0;
    prepareData_uint8(&cur_G, cur_A, A, G_len, A_len, cur_test, "", 0);

    union gen cur_gen;
    for (int i = 0; i < 8; i++){
        cur_gen.bin = cur_G[i];
    }
     */

    // Test mapPhenotypeValuesToShort
    /*
    char *phn_name = malloc(25 * sizeof(*phn_name));
    phn_name[0] = 'w'; phn_name[1] = 'o'; phn_name[2] = 'r'; phn_name[3] = 'd'; phn_name[4] = '1';
    phn_name[5] = 'w'; phn_name[6] = 'o'; phn_name[7] = 'r'; phn_name[8] = 'd'; phn_name[9] = '2';
    phn_name[10] = 'w'; phn_name[11] = 'o'; phn_name[12] = 'r'; phn_name[13] = 'd'; phn_name[14] = '1';
    phn_name[15] = 'w'; phn_name[16] = 'o'; phn_name[17] = 'r'; phn_name[18] = 'd'; phn_name[19] = '4';
    phn_name[20] = 'w'; phn_name[21] = 'o'; phn_name[22] = 'r'; phn_name[23] = 'd'; phn_name[24] = '2';
    size_t phn_cnt = 5;
    ptrdiff_t phn_off = 5;
    unsigned short * new_phn = mapPhenotypeValuesToShort(phn_name, phn_cnt, phn_off);
    */


    // Test calcPValue_n
    /*
    size_t phn_len = 6;

    uint8_t *genotype = malloc(sizeof(uint8_t) * 9);
    union gen gen;
    gen.n0 = 0; gen.n1 = 1; gen.n2 = 0; gen.n3 = 0;
    genotype[0] = gen.bin;
    gen.n0 = 0; gen.n1 = 0; gen.n2 = 0; gen.n3 = 1;
    genotype[1] = gen.bin;
    gen.n0 = 1; gen.n1 = 1; gen.n2 = 0; gen.n3 = 0;
    genotype[2] = gen.bin;

    gen.n0 = 0; gen.n1 = 1; gen.n2 = 0; gen.n3 = 0;
    genotype[3] = gen.bin;
    gen.n0 = 1; gen.n1 = 0; gen.n2 = 0; gen.n3 = 1;
    genotype[4] = gen.bin;
    gen.n0 = 1; gen.n1 = 0; gen.n2 = 1; gen.n3 = 1;
    genotype[5] = gen.bin;

    gen.n0 = 0; gen.n1 = 0; gen.n2 = 0; gen.n3 = 1;
    genotype[6] = gen.bin;
    gen.n0 = 0; gen.n1 = 0; gen.n2 = 0; gen.n3 = 1;
    genotype[7] = gen.bin;
    gen.n0 = 1; gen.n1 = 1; gen.n2 = 1; gen.n3 = 1;
    genotype[8] = gen.bin;


    //bool res = checkNumElemInGenotype2(genotype, phn_len);
    //if (res){
        //printf("Ok!");
    //}

    unsigned short *phenotype = malloc(sizeof(*phenotype) * 12);
    phenotype[0] = 1; phenotype[1] = 3; phenotype[2] = 0; phenotype[3] = 2; phenotype[4] = 1; phenotype[5] = 0;
    phenotype[6] = 1; phenotype[7] = 3; phenotype[8] = 0; phenotype[9] = 2; phenotype[10] = 1; phenotype[11] = 0;


    size_t G_len = 36, A_len = 12;

    double pvalres = calcPValue_n(genotype, phenotype,  G_len,  A_len);

    free(genotype);
    free(phenotype);
    */

    // Test the whole thing
    TableEntry *tests = (TableEntry *) malloc (sizeof(TableEntry)); // Just one entry for now
    //tests[0] = {.test_index = "1", .chromosome_index = "1", .ID = "cd", .lower = 0, .upper = 3, .adjusted_p_value = 0.0};
    tests[0].test_index = "1";
    tests[0].chromosome_index = "1";
    tests[0].ID = "cd"; // r, d, a, d results: 0.7 1.0 0.6 0.2 0.0
    tests[0].lower = 0;
    tests[0].upper = 3;
    tests[0].adjusted_p_value = 0.0;
    size_t tests_len = 1;

    unsigned short *phen = (unsigned short *) malloc (6 * sizeof(unsigned short));
    phen[0] = 1; phen[1] = 3; phen[2] = 0; phen[3] = 2; phen[4] = 1; phen[5] = 0;
    size_t phen_len = 6;

    ExecutionParameters parameters;
    parameters.maxReplications = 10; // Should be around 10^(-8)
    parameters.k = 10;
    parameters.isAdaptive = 0;

    char *path = "does/not/matter/now";

    int type = 0;

    adjustPValue(tests, phen, tests_len, phen_len, parameters, path, type);// tests, phen, tests_len, phen_len, parameters, path, type

    for (int i = 0; i < tests_len; ++i) {
        printf("%f\n", tests[i].adjusted_p_value);
    }

    return 0;
}