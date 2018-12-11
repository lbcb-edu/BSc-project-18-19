#include "gtest/gtest.h"
#include "white_alignment.h"
#include "white_minimizers.h"

namespace {

	TEST (WhiteAlignmentTest, GlobalAlignmentWorks) {
		std::string q = "GGGCCCATAGCTCAGTGGTAGAGTGCCTCCTTTGCAAGGAGGATGCCCTGGGTTCGAATCCCAGTGGGTCCAT";
		std::string t = "GGGCCCATAGCTCAGCCTGGGAGAGCGCCGCCCTTGCAAGGCGGAGGCCCCGGGTTCAAATCCCGGTGGGTCCAT";

		int match = 2;
		int mismatch = -2;
		int gap = -2;

		EXPECT_EQ (106, white::pairwise_alignment (q.c_str(), q.size(), t.c_str(), t.size(), white::AlignmentType::global, match, mismatch, gap));
	}

	TEST (WhiteAlignmentTest, SemiGlobalAlignmentWorks) {
		std::string q = {"ATTGGAA"};
        	std::string t = {"CCCACTTTTTGG"};

        	int match = 4;
        	int mismatch = -1;
        	int gap = -2;

		std::string cigar;
  		unsigned int target_begin = 0;

  		int value = white::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), white::AlignmentType::semi_global, match, mismatch, gap, cigar, target_begin);

		EXPECT_EQ(cigar, "1X4=2S");
		EXPECT_EQ(target_begin, 8);
	}

	TEST (WhiteAlignmentTest, LocalAlignmentWorks) {
        	std::string s1 = "CCAATGCCACAAAACATCTGTCTCTAACTGGTGTGTGTGT";
        	std::string s2 = "CCAGCCCAAAATCTGTTTTAATGGTGGATTTGTGT";

        	int match = 2;
        	int mismatch = -2;
      		int gap = -2;

        	EXPECT_EQ (40, white::pairwise_alignment (s1.c_str(), s1.size(), s2.c_str(), s2.size(), white::AlignmentType::local, match, mismatch, gap));
	}

	TEST (WhiteMinimizersTest, DeterminingMinimizersWorks) {
	std::string s1 = "TCAGGAAGAAGCAGA";

	int k = 3;
	int window_length = 4;

	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;

	minimizers.push_back(std::make_tuple(133, 2, true));
	minimizers.push_back(std::make_tuple(113, 5, true));
	minimizers.push_back(std::make_tuple(113, 8, true));
	minimizers.push_back(std::make_tuple(131, 12, true));

	EXPECT_EQ (minimizers, white::minimizers(s1.c_str(), s1.size(), k, window_length));

	std::string s2 = "GTCATGCACGTTCAC";

        k = 3;
        window_length = 4;

	minimizers.clear();
        minimizers.push_back(std::make_tuple(112, 3, false));
        minimizers.push_back(std::make_tuple(123, 7, true));
        minimizers.push_back(std::make_tuple(143, 10, false));

        EXPECT_EQ (minimizers, white::minimizers(s2.c_str(), s2.size(), k, window_length));

	std::string s3 = "GTCATGCACGTTCAC";

        k = 6;
        window_length = 5;

	minimizers.clear();
        minimizers.push_back(std::make_tuple(112343, 3, false));
        minimizers.push_back(std::make_tuple(123432, 4, false));
	minimizers.push_back(std::make_tuple(123442, 7, true));

        EXPECT_EQ (minimizers, white::minimizers(s3.c_str(), s3.size(), k, window_length));
	}

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
