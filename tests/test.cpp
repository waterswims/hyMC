const double pi = 3.141592653589793;

#include "mklrand_test.hpp"
#include "1d_hamil_tests.hpp"
#include "./leapfrog_test.hpp"
#include "./hmc_test.hpp"
#include "gtest/gtest.h"

// Run all tests
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
