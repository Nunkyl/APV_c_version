//
// Created by doby on 17/11/17.
//

#include "PValue.h"
#include "stdlib.h"
#include <stdio.h>
#include <gsl/gsl_sf_gamma.h>
#include <string.h>
#include <time.h>
#include <stddef.h>


void adjustPValue(TableEntry *tests, const unsigned short *A, size_t tests_len, size_t A_len, const ExecutionParameters cont, char *path, int type) { // TableEntry *tests, const unsigned short *A, size_t tests_len, size_t A_len, const ExecutionParameters cont, char *path, int type
    //srand ((int) time(NULL));

    uint8_t *cur_G;
    unsigned short *cur_A;
    double D_main = 0.0, D_cur = 0.0;
    int s, m, k, all_iter;
    size_t cur_G_len, cur_A_len, G_len = (A_len+3)/4 * (tests->upper - tests->lower + 1); // G_len - size in bytes
    if (tests->ID == "a"){
        cur_G_len = 2 * G_len;
        cur_A_len = 2 * A_len;
    }
    else {
        cur_G_len = G_len;
        cur_A_len = A_len;
    }

    if (cont.isAdaptive) k = cont.k;
    else k = cont.maxReplications;

    //#pragma omp parallel for
    for (size_t i = 0; i < tests_len; ++i) {
        cur_A_len = A_len;
        cur_G = malloc (2 * G_len * sizeof(*cur_G));
        cur_A = malloc (2 * A_len * sizeof(*cur_A));
        prepareData(&cur_G, cur_A, A, cur_G_len, cur_A_len, tests[i], path, type);
        D_main = calcPValue(cur_G, cur_A, cur_G_len, cur_A_len);
        s = 0; m = 0; all_iter = 0;
        while (s < cont.maxReplications && m < k && all_iter < cont.maxReplications){
            all_iter++;
            phenotypeRandomPermutation(cur_A, cur_A_len);
            D_cur = calcPValue(cur_G, cur_A, cur_G_len, cur_A_len);
            if (isnan(D_cur)) continue; // If D_cur is NaN go to the next iteration ((D_cur != D_cur))
            s++;
            if (D_cur > D_main) m++;
        }
        if (s == 0) tests[i].adjusted_p_value = -1;
        else tests[i].adjusted_p_value = (double)m/(double)s;
        free(cur_G);
        free(cur_A);
    }

}

void prepareData(uint8_t **cur_G, unsigned short *cur_A, const unsigned short *A, size_t G_len, size_t A_len, TableEntry cur_test, char *path, int type){

    union gen cur_gen;
    unsigned short bits[4];

    size_t row_num = cur_test.upper - cur_test.lower + 1;
    createGenotypeMatrix(*cur_G, path, G_len, cur_test.lower, cur_test.upper, type); // Read the appropriate part of the genotype

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

double calcPValue(const uint8_t *cur_G, const unsigned short *cur_A, size_t G_len, size_t A_len){

    double D = 0.0;     // Return value
    int D_num = 0;      // Amount of valid D values
    int V_rows, V_cols; // Size of V matrix
    int G_car;          // The cardinality of the sets of values of the phenotype and genotype
    int s;              // Parameters for calculating the gamma function
    double chi_sqr = 0.0;
    size_t row_num = G_len/((A_len+3)/4); // The amount of rows in the current genotype matrix
    PhenotypeStatistics phen_stat;
    size_t row_len_bytes = (A_len + 3)/4;

    size_t snp_len = (A_len + 3)/4;// Length of line in genotype matrix

    phen_stat = calcNumElem(cur_A, A_len);
    if (phen_stat.num_elements <= 1) return NAN; //?????????? quite NaN

    //Allocate memory for V matrix here, it will always be 4x(phen_stat.max_element+1)
    int *V = malloc(4*(phen_stat.max_element+1)*sizeof(*V));
    if (!V) { /* Error handling */  }

    for (int i = 0; i < row_num; ++i) {

        if (!checkNumElemInGenotype(cur_G + i*(A_len+3)/4, A_len)) continue;
        fillVMatrix(cur_G + i*row_len_bytes,  cur_A, V, 4, (phen_stat.max_element+1), A_len);
        G_car = calcNumElementsInGenotype(V, 4, (phen_stat.max_element+1)); // Check this!

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

void createGenotypeMatrix(uint8_t *genotype, char *path, size_t G_len, size_t lower, size_t upper, int type){
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


        /*
        gen[0]  = 0; gen[1]  = 0; gen[2]  = 0; gen[3]  = 1; gen[4]  = 1; gen[5]  = 0;
        gen[6]  = 2; gen[7]  = 1; gen[8]  = 0; gen[9]  = 0; gen[10] = 2; gen[11] = 1;
        gen[12] = 1; gen[13] = 0; gen[14] = 0; gen[15] = 2; gen[16] = 2; gen[17] = 1;
        gen[18] = 0; gen[19] = 0; gen[20] = 1; gen[21] = 1; gen[22] = 1; gen[23] = 0;
        */

        /*
        gen[0]  = 2; gen[1]  = 2; gen[2]  = 1; gen[3]  = 0; gen[4]  = 1; gen[5]  = 2;
        gen[6]  = 0; gen[7]  = 1; gen[8]  = 0; gen[9]  = 2; gen[10] = 2; gen[11] = 1;
        gen[12] = 2; gen[13] = 0; gen[14] = 1; gen[15] = 0; gen[16] = 1; gen[17] = 2;
        */

        /*
        gen[0]  = 0; gen[1]  = 3; gen[2]  = 0; gen[3]  = 1; gen[4]  = 1; gen[5]  = 2;
        gen[6]  = 1; gen[7]  = 0; gen[8]  = 2; gen[9]  = 2; gen[10] = 1; gen[11] = 1;
         */


        /*
        gen[0]  = 0; gen[1]  = 2; gen[2]  = 1; gen[3]  = 1; gen[4]  = 0; gen[5]  = 0;
        gen[6]  = 0; gen[7]  = 2; gen[8]  = 1; gen[9]  = 0; gen[10] = 2; gen[11] = 1;
        gen[12] = 0; gen[13] = 1; gen[14] = 1; gen[15] = 2; gen[16] = 1; gen[17] = 1;
        */

        /*
        gen[0]  = 2; gen[1]  = 0; gen[2]  = 0; gen[3]  = 0; gen[4]  = 2; gen[5]  = 1;
        gen[6]  = 3; gen[7]  = 0; gen[8]  = 3; gen[9]  = 1; gen[10] = 1; gen[11] = 2;
        gen[12] = 0; gen[13] = 1; gen[14] = 2; gen[15] = 0; gen[16] = 2; gen[17] = 1;
        gen[18] = 0; gen[19] = 1; gen[20] = 3; gen[21] = 2; gen[22] = 0; gen[23] = 0;
         */


    }
    if (type == 1) {} // If data comes from text file
    if (type == 2) {} // If data comes from PLINK
    //return gen;
}

bool checkNumElemInGenotype(const uint8_t *genotype, size_t phn_len) {
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
        if(j == (phn_len + 3)/4 - 1  && phn_len%4 != 0) bits_num = phn_len%4;
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

bool checkNumElemInGenotype2(const uint8_t *genotype, size_t phn_len) {
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
        if(j == (phn_len + 3)/4 - 1) bits_num = phn_len%4 + 1;
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

AlternativeHypothesisType hashIt (const char *inString) {
    if (inString == "cd") return eCD;
    if (inString == "r")  return eR;
    if (inString == "d")  return eD;
    if (inString == "a")  return eA;
    return eCD;
}

PhenotypeStatistics calcNumElem(const unsigned short *phenotype, size_t phen_len){
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

/*
unsigned short * doubleSizeOfPhenotype(unsigned short *cur_A, size_t A_len){
    size_t size_cur_A =  A_len * sizeof(*cur_A); // The amount of memory required for the old phenotype array
    unsigned short *buf_A =  malloc(2 * size_cur_A);
    if (!buf_A) { }
    memcpy(buf_A, cur_A, size_cur_A); // Copy cur_A into the first half of the new array
    memcpy(buf_A + A_len, cur_A, size_cur_A); // Copy cur_A into the second half of the new array
    free(cur_A); // ????
    return buf_A;
}


unsigned short * doubleSizeOfGenotype(unsigned short *cur_G, size_t G_len, size_t row_num){
    size_t col_num = G_len/row_num;
    size_t size_cur_G = G_len * sizeof(*cur_G);
    size_t row_len = col_num * sizeof(*cur_G);
    unsigned short *buf_G =  malloc(2 * size_cur_G);
    if (!buf_G) {  }

    for (size_t i = 0; i < row_num; i++){
        memcpy(buf_G + col_num*i*2, cur_G + col_num*i, row_len);
        memcpy(buf_G + col_num*(i*2+1), cur_G + col_num*i, row_len);
    }
    free(cur_G); // ????
    return buf_G;
}
*/


void fillVMatrix(const uint8_t *cur_G, const unsigned short *cur_A, int *V, int V_rows, int V_cols, size_t col_num){  // col_num - The amount of patients

    union gen cur_gen;
    uint8_t bits[4];
    int bits_num = 4;

    // Fill V with zeros
    memset (V,0,V_rows * V_cols * sizeof(*V));

    // Fill V (G_car x A_car)
    for (size_t j = 0; j < (col_num + 3)/4; ++j) {
        cur_gen.bin = cur_G[j]; // uint8_t a = cur_G[j];
        bits[0] = cur_gen.n0, bits[1] = cur_gen.n1, bits[2] = cur_gen.n2, bits[3] = cur_gen.n3;
        if(j == (col_num + 3)/4 - 1 && col_num%4 != 0) bits_num = col_num%4;
        for (int i = 0; i < bits_num; i++) {
            V[bits[i] * V_cols + cur_A[j*4 + i]]++; // bits => a & 3; a >>= 2
        }
    }
}

double calculateChiSqr(const int *V, int V_rows, int V_cols){

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

size_t rand_interval(size_t min, size_t max)
{
    size_t r;
    size_t range = 1 + max - min;
    size_t buckets = RAND_MAX / range;
    size_t limit = buckets * range;

    //srand (time(NULL));

    do {
        r = (size_t)rand();
    } while (r >= limit);

    return min + (size_t)(r / buckets);
}

// Permutate phenotype
void phenotypeRandomPermutation(unsigned short *phenotype, size_t phenotypeLength) {
    for (unsigned i = 0; i < phenotypeLength; i++) { //  Knuth shuffle
        size_t k = rand_interval(i, phenotypeLength-1);
        unsigned short tmp = phenotype[k];
        phenotype[k] = phenotype[i];
        phenotype[i] = tmp;
    }
}

// Calculates the exact number of elements in a genotype vector based on the V matrix
int calcNumElementsInGenotype(const int *V, int rowNum, int colNum){
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


unsigned short * mapPhenotypeValuesToShort(const char *phn_name, size_t phn_cnt, ptrdiff_t phn_off){
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