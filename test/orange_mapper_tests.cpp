#include "gtest/gtest.h"
#include "orange_alignment.h"

TEST (Pairwise_alignment, global) {
	
	std::string cigar;
	unsigned int target_start = 0;

	EXPECT_EQ (orange::pairwise_alignment ("GCATGCU", 7, "GATTACA", 7, orange::AlignmentType::global, 1, -1, -1, cigar, target_start), 0);
	EXPECT_EQ (target_start, 0);
}

TEST (Pairwise_alignment, local) {
	std::string cigar;
	unsigned int target_start = 0;

	EXPECT_EQ (orange::pairwise_alignment ("ACCTAAGG", 8, "GGCTCAATCA", 10, orange::AlignmentType::local, 2, -1, -2, cigar, target_start), 6);
	EXPECT_EQ (target_start, 2);
}

TEST (Pairwise_alignment, semi_global) {
	std::string cigar;
	unsigned int target_start = 0;

	EXPECT_EQ (orange::pairwise_alignment ("ACCCAAGGG", 9, "GGCTCCATTA", 10, orange::AlignmentType::semi_global, 1, -1, -1, cigar, target_start), 1);
	EXPECT_EQ (target_start, 1);
}
