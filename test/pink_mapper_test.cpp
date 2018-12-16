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

    EXPECT_EQ(std::get<0>(minimizers.at(0)), 14);
    EXPECT_EQ(std::get<1>(minimizers.at(0)), 2);
    EXPECT_EQ(std::get<2>(minimizers.at(0)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(1)), 14);
    EXPECT_EQ(std::get<1>(minimizers.at(1)), 3);
    EXPECT_EQ(std::get<2>(minimizers.at(1)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(2)), 17);
    EXPECT_EQ(std::get<1>(minimizers.at(2)), 6);
    EXPECT_EQ(std::get<2>(minimizers.at(2)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(3)), 6);
    EXPECT_EQ(std::get<1>(minimizers.at(3)), 7);
    EXPECT_EQ(std::get<2>(minimizers.at(3)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(4)), 1);
    EXPECT_EQ(std::get<1>(minimizers.at(4)), 9);
    EXPECT_EQ(std::get<2>(minimizers.at(4)), 1);
    EXPECT_EQ(std::get<0>(minimizers.at(5)), 17);
    EXPECT_EQ(std::get<1>(minimizers.at(5)), 12);
    EXPECT_EQ(std::get<2>(minimizers.at(5)), 0);
    EXPECT_EQ(std::get<0>(minimizers.at(6)), 33);
    EXPECT_EQ(std::get<1>(minimizers.at(6)), 0);
    EXPECT_EQ(std::get<2>(minimizers.at(6)), 1);
}
