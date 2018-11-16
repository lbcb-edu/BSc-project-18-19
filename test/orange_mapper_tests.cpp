#include "gtest/gtest.h"
#include "orange_alignment.h"

TEST (Pairwise_alignment, global) {
    EXPECT_EQ (orange::pairwise_alignment ("GCATGCU", 7, "GATTACA", 7, orange::AlignmentType::global, 1, -1, -1), 0);
}

TEST (Pairwise_alignment, local) {
    EXPECT_EQ (orange::pairwise_alignment ("ACCTAAGG", 8, "GGCTCAATCA", 10, orange::AlignmentType::local, 2, -1, -2), 6);
}