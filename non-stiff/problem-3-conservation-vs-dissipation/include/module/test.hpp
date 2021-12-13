#pragma once
#include "function/prec.hpp"

static int test_err = 0;
static int num_test_pass = 0;
static int num_test_fail = 0;

#define test(cond) test_((cond), #cond, __LINE__)
void test_(const bool cond, const char *label, const int line) {
        const char *stat;
        if (cond) {
                stat = "Passed";
                num_test_pass += 1;
        } else {
                test_err = 1;
                stat = "Failed";
                num_test_fail += 1;
        }

        printf("\t %d:%s \t %s \n", line, label, stat);
}

void test_stats() {
        printf("%d out of %d test(s) passed\n", num_test_pass,
               num_test_fail + num_test_pass);
}
