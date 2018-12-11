#include "gtest/gtest.h"
#include <pink_alignment.hpp>
#include <pink_minimizers.hpp>

TEST (Pairwise_alignment, global) {
    EXPECT_EQ (pink::pairwise_alignment ("TGCATAT", 7, "ATCCGAT", 7, pink::global, 1, -1, -1), 0);
}

TEST (Pairwise_alignment, semi_global) {
    EXPECT_EQ (pink::pairwise_alignment ("TGCATAT", 7, "ATCCGAT", 7, pink::semi_global, 1, -1, -1), 2);
}

TEST (Pairwise_alignment, local) {
    EXPECT_EQ (pink::pairwise_alignment ("TGCATAT", 7, "ATCCGAT", 7, pink::local, 1, -1, -1), 2);
}

TEST (Minimizers, TestAGTCCTTCTTCGTCT) {
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = pink::minimizers("AGTCCTTCTTCGTCT", 15, 3, 3);

    EXPECT_EQ(std::get<0>(minimizers.at(0)), 32);
    EXPECT_EQ(std::get<1>(minimizers.at(0)), 2);
    EXPECT_EQ(std::get<2>(minimizers.at(0)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(1)), 32);
    EXPECT_EQ(std::get<1>(minimizers.at(1)), 2);
    EXPECT_EQ(std::get<2>(minimizers.at(1)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(2)), 32);
    EXPECT_EQ(std::get<1>(minimizers.at(2)), 2);
    EXPECT_EQ(std::get<2>(minimizers.at(2)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(3)), 32);
    EXPECT_EQ(std::get<1>(minimizers.at(3)), 3);
    EXPECT_EQ(std::get<2>(minimizers.at(3)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(4)), 101);
    EXPECT_EQ(std::get<1>(minimizers.at(4)), 6);
    EXPECT_EQ(std::get<2>(minimizers.at(4)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(5)), 12);
    EXPECT_EQ(std::get<1>(minimizers.at(5)), 7);
    EXPECT_EQ(std::get<2>(minimizers.at(5)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(6)), 12);
    EXPECT_EQ(std::get<1>(minimizers.at(6)), 7);
    EXPECT_EQ(std::get<2>(minimizers.at(6)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(7)), 1);
    EXPECT_EQ(std::get<1>(minimizers.at(7)), 9);
    EXPECT_EQ(std::get<2>(minimizers.at(7)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(8)), 1);
    EXPECT_EQ(std::get<1>(minimizers.at(8)), 9);
    EXPECT_EQ(std::get<2>(minimizers.at(8)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(9)), 1);
    EXPECT_EQ(std::get<1>(minimizers.at(9)), 9);
    EXPECT_EQ(std::get<2>(minimizers.at(9)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(10)), 101);
    EXPECT_EQ(std::get<1>(minimizers.at(10)), 12);
    EXPECT_EQ(std::get<2>(minimizers.at(10)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(11)), 231);
    EXPECT_EQ(std::get<1>(minimizers.at(11)), 0);
    EXPECT_EQ(std::get<2>(minimizers.at(11)), 0);
}

TEST(MapperTest, TestCGACT){
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = pink::minimizers("CGACT", 5, 3, 1);
    
    EXPECT_EQ(std::get<0>(minimizers.at(0)), 100);
    EXPECT_EQ(std::get<1>(minimizers.at(0)), 0);
    EXPECT_EQ(std::get<2>(minimizers.at(0)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(1)), 310);
    EXPECT_EQ(std::get<1>(minimizers.at(1)), 1);
    EXPECT_EQ(std::get<2>(minimizers.at(1)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(2)), 201);
    EXPECT_EQ(std::get<1>(minimizers.at(2)), 2);
    EXPECT_EQ(std::get<2>(minimizers.at(2)), 0);
}
