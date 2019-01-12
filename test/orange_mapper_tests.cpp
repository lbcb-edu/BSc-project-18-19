#include "gtest/gtest.h"
#include "orange_alignment.h"
#include "orange_minimizers.h"

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

TEST (Minimizers, 3mer_windowLenght3) {
	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizer_vec;

	minimizer_vec.emplace_back(18, 0, 1);
	minimizer_vec.emplace_back(11, 1, 1);
	minimizer_vec.emplace_back(6, 4, 1);
	minimizer_vec.emplace_back(6, 7, 0);
	minimizer_vec.emplace_back(2, 10, 1);
	minimizer_vec.emplace_back(11, 11, 1);
	minimizer_vec.emplace_back(17, 12, 0);

	std::vector<std::tuple<unsigned int, unsigned int, bool>> res_vec = orange::minimizers("TGACGTACATGGACA", 15, 3, 3);
	std::sort(res_vec.begin(), res_vec.end(), [](const std::tuple<unsigned int, unsigned int, bool> &a, const std::tuple<unsigned int, unsigned int, bool> &b)
							{return std::get<1>(a) < std::get<1>(b);});

	EXPECT_EQ (res_vec, minimizer_vec);
}
