#include "blue_alignment.hpp"
#include "blue_minimizers.hpp"
#include "gtest/gtest.h"
#include <string>
#include <vector>
#include <tuple>

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

        EXPECT_EQ(4, blue::pairwise_alignment(query, 4, target, 8, blue::semi_global, 1, -2, -1, cigar, target_begin));
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

    TEST(BlueMinimizersTest, VectorLengthTest) {
        char * sequence = "TGACGTACAT";

        std::vector<std::tuple<unsigned int, unsigned int, bool>> result = blue::minimizers(sequence, 10, 3, 4);

        EXPECT_EQ(6, result.size());

    }

    TEST(BlueMinimizersTest, VectorElementTest) {
        char * sequence = "TAACGTG";

        std::vector<std::tuple<unsigned int, unsigned int, bool>> result = blue::minimizers(sequence, 7, 3, 3);

        EXPECT_EQ(std::make_tuple(4, 2, 1), result.at(0));

    }

}

