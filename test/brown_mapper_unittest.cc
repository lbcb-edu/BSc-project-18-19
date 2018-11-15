#include "gtest/gtest.h"
#include <string>
#include "brown_alignment.hpp"

namespace {

	TEST(MapperTest, SemiGlobalTest){
		std::string q = {"ATTGGAA"};
  		std::string t = {"CCCACTTTTTGG"};
		std::string cigar;
  		unsigned int target_begin = 0;
  		int value = brown::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), brown::AlignmentType::semi_global, 4, -1, -2, cigar, target_begin);

		EXPECT_EQ(cigar, "1X4=2S");
		EXPECT_EQ(target_begin, 7);
	}

	TEST(MapperTest, LocalTest){
		std::string q = {"TTCCGCCAA"};
  		std::string t = {"AACCCCTT"};
  		std::string cigar;
  		unsigned int target_begin = 0;
  		int value = brown::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), brown::AlignmentType::local, 4, -1, -2, cigar, target_begin);

		EXPECT_EQ(cigar, "2S2=1I2=2S");
		EXPECT_EQ(target_begin, 2);
	}

	TEST(MapperTest, GlobalTest){
		std::string q = {"TGCATAT"};
  		std::string t = {"ATCCGAT"};
  		std::string cigar;
  		unsigned int target_begin = 0;
  		int value = brown::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), brown::AlignmentType::global, 4, -1, -2, cigar, target_begin);

		EXPECT_EQ(cigar, "1D1=1I1=2X2=");
		EXPECT_EQ(target_begin, 0);
	}
}

// int main(int argc, char **argv) {
//   ::testing::InitGoogleTest(&argc, argv);
//   return RUN_ALL_TESTS();
// }