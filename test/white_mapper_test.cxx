#include "gtest/gtest.h"
#include "white_alignment.h"

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

		EXPECT_EQ(cigar, "1X4M2S");
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

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
