#include "blue_alignment.hpp"
#include "gtest/gtest.h"
#include <string>

namespace {
    TEST(BlueAlignmentTest, GlobalAlignmentTest){
        char * query = "CTCTGTTCG";
        char * target = "CGTATCTTGA";

        EXPECT_EQ(-5, blue::pairwise_alignment(query, 9, target, 10, blue::global, 0, -1, -1));
    }

    TEST(BlueAlignmentTest, SemiGlobalAlignmentTest){
        char * query = "TCCG";
        char * target = "ACTCCGAT";

        std::string cigar;
        unsigned int target_begin = 0;

        EXPECT_EQ(2, blue::pairwise_alignment(query, 9, target, 10, blue::semi_global, 1, -2, -1, cigar, target_begin));
        EXPECT_EQ(target_begin, 2);
    }

    TEST(BlueAlignmentTest, LocalAlignmentTest) {
        char * query = "ACCTAAGG";
        char * target = "GGCTCAATCA";

        std::string cigar;
        unsigned int target_begin = 0;

        EXPECT_EQ(6, blue::pairwise_alignment(query, 8, target, 10, blue::local, 2, -1, -2, cigar, target_begin));
        EXPECT_EQ(target_begin, 2);
    }

}

