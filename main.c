#include <stdio.h>
#include <stdlib.h>
#include "PValue.h"
#include <gsl/gsl_sf_gamma.h>

int main() {

    // Run unit tests
    //testing::InitGoogleTest(&argc, argv);
    //RUN_ALL_TESTS();


    TableEntry *tests = (TableEntry *) malloc (sizeof(TableEntry)); // Just one entry for now
    //tests[0] = {.test_index = "1", .chromosome_index = "1", .ID = "cd", .lower = 0, .upper = 3, .adjusted_p_value = 0.0};
    tests[0].test_index = "1";
    tests[0].chromosome_index = "1";
    tests[0].ID = "cd";
    tests[0].lower = 0;
    tests[0].upper = 2;
    tests[0].adjusted_p_value = 0.0;
    size_t tests_len = 1;

    unsigned short *phen = (unsigned short *) malloc (6 * sizeof(unsigned short));
    phen[0] = 1; phen[1] = 3; phen[2] = 0; phen[3] = 2; phen[4] = 1; phen[5] = 0;
    size_t phen_len = 6;

    ExecutionParameters parameters;
    parameters.maxReplications = 10; // Should be around 10^(-8)
    parameters.k = 10;
    parameters.isAdaptive = 0;

    char *path = "does/not/matter/now";//

    int type = 0;

    adjustPValue(tests, phen, tests_len, phen_len, parameters, path, type);// tests, phen, tests_len, phen_len, parameters, path, type

    //void adjustPValue(TableEntry *tests, const unsigned short *A, size_t tests_len, size_t A_len, const ExecutionParameters cont, char *path, int type)

    for (int i = 0; i < tests_len; ++i) {
        printf("%f\n", tests[i].adjusted_p_value);
    }

    return 0;
}