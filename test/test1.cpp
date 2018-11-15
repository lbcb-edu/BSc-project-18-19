#include <src/pink_alignment.hpp>

TEST (Pairwise_alignment, semi_global) {
    EXPECT_EQ (2, pairwise_alignment ("TGCATAT", 7, "ATCCGAT", 7, semi_global, 1, -1, -1));
}
TEST (Pairwise_alignment, global) {
    EXPECT_EQ (0, pairwise_alignment ("TGCATAT", 7, "ATCCGAT", 7, global, 1, -1, -1));
}
TEST (Pairwise_alignment, local) {
    EXPECT_EQ (2, pairwise_alignment ("TGCATAT", 7, "ATCCGAT", 7, local, 1, -1, -1));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
