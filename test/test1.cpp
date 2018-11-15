#include "gtest/gtest.h"
#include <src/pink_alignment.hpp>


TEST (Pairwise_alignment, semi_global) {
    EXPECT_EQ (pink::pairwise_alignment ("TGCATAT", 7, "ATCCGAT", 7, semi_global, 1, -1, -1), 2);
}
TEST (Pairwise_alignment, global) {
    EXPECT_EQ (pink::pairwise_alignment ("TGCATAT", 7, "ATCCGAT", 7, global, 1, -1, -1), 0);
}
TEST (Pairwise_alignment, local) {
    EXPECT_EQ (pink::pairwise_alignment ("TGCATAT", 7, "ATCCGAT", 7, local, 1, -1, -1), 2);
}

//int main(int argc, char **argv) {
//    ::testing::InitGoogleTest(&argc, argv);
//    return RUN_ALL_TESTS();
//}
