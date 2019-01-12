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

TEST (Minimizers, TestTGACGTACATGGACA) {
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = pink::minimizers("TGACGTACATGGACA", 15, 3, 3);

    EXPECT_EQ(std::get<0>(minimizers.at(0)), 17);
    EXPECT_EQ(std::get<1>(minimizers.at(0)), 12);
    EXPECT_EQ(std::get<2>(minimizers.at(0)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(1)), 33);
    EXPECT_EQ(std::get<1>(minimizers.at(1)), 0);
    EXPECT_EQ(std::get<2>(minimizers.at(1)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(2)), 1);
    EXPECT_EQ(std::get<1>(minimizers.at(2)), 9);
    EXPECT_EQ(std::get<2>(minimizers.at(2)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(3)), 6);
    EXPECT_EQ(std::get<1>(minimizers.at(3)), 8);
    EXPECT_EQ(std::get<2>(minimizers.at(3)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(4)), 6);
    EXPECT_EQ(std::get<1>(minimizers.at(4)), 7);
    EXPECT_EQ(std::get<2>(minimizers.at(4)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(5)), 17);
    EXPECT_EQ(std::get<1>(minimizers.at(5)), 6);
    EXPECT_EQ(std::get<2>(minimizers.at(5)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(6)), 14);
    EXPECT_EQ(std::get<1>(minimizers.at(6)), 3);
    EXPECT_EQ(std::get<2>(minimizers.at(6)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(7)), 14);
    EXPECT_EQ(std::get<1>(minimizers.at(7)), 2);
    EXPECT_EQ(std::get<2>(minimizers.at(7)), 1);
    
}
