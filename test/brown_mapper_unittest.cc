#include "gtest/gtest.h"
#include <string>
#include "brown_alignment.hpp"
#include "brown_minimizers.hpp"

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

		TEST(MapperTest, TestCGAC){
 		std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = brown::minimizers("CGAC",4, 3 , 1);

		EXPECT_EQ(std::get<0>(minimizers.at(0)), 13);
		EXPECT_EQ(std::get<1>(minimizers.at(0)), 0);
		EXPECT_EQ(std::get<2>(minimizers.at(0)), 0);
		EXPECT_EQ(std::get<0>(minimizers.at(1)), 52);
		EXPECT_EQ(std::get<1>(minimizers.at(1)), 1);
		EXPECT_EQ(std::get<2>(minimizers.at(1)), 0);
		EXPECT_EQ(minimizers.size(), 2);
	}

	TEST(MapperTest, TestCGACT){
		std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = brown::minimizers("CGACT", 5, 3, 1);

		EXPECT_EQ(std::get<0>(minimizers.at(0)), 13);
		EXPECT_EQ(std::get<1>(minimizers.at(0)), 0);
		EXPECT_EQ(std::get<2>(minimizers.at(0)), 0);
		EXPECT_EQ(std::get<0>(minimizers.at(1)), 52);
		EXPECT_EQ(std::get<1>(minimizers.at(1)), 1);
		EXPECT_EQ(std::get<2>(minimizers.at(1)), 0);
		EXPECT_EQ(std::get<0>(minimizers.at(2)), 18);
		EXPECT_EQ(std::get<1>(minimizers.at(2)), 2);
		EXPECT_EQ(std::get<2>(minimizers.at(2)), 0);
		EXPECT_EQ(minimizers.size(), 3);
	}

	TEST(MapperTest, TestTGACGTACATGGACA){
		std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = brown::minimizers("TGACGTACATGGACA", 15, 3, 3);

		EXPECT_EQ(std::get<0>(minimizers.at(0)), 33);
		EXPECT_EQ(std::get<1>(minimizers.at(0)), 0);
		EXPECT_EQ(std::get<2>(minimizers.at(0)), 1);
		EXPECT_EQ(std::get<0>(minimizers.at(1)), 14);
		EXPECT_EQ(std::get<1>(minimizers.at(1)), 2);
		EXPECT_EQ(std::get<2>(minimizers.at(1)), 1);
		EXPECT_EQ(std::get<0>(minimizers.at(2)), 14);
		EXPECT_EQ(std::get<1>(minimizers.at(2)), 3);
		EXPECT_EQ(std::get<2>(minimizers.at(2)), 0);
		EXPECT_EQ(std::get<0>(minimizers.at(3)), 17);
		EXPECT_EQ(std::get<1>(minimizers.at(3)), 6);
		EXPECT_EQ(std::get<2>(minimizers.at(3)), 0);
		EXPECT_EQ(std::get<0>(minimizers.at(4)), 6);
		EXPECT_EQ(std::get<1>(minimizers.at(4)), 7);
		EXPECT_EQ(std::get<2>(minimizers.at(4)), 0);
		EXPECT_EQ(std::get<0>(minimizers.at(5)), 6);
		EXPECT_EQ(std::get<1>(minimizers.at(5)), 8);
		EXPECT_EQ(std::get<2>(minimizers.at(5)), 1);
		EXPECT_EQ(std::get<0>(minimizers.at(6)), 1);
		EXPECT_EQ(std::get<1>(minimizers.at(6)), 9);
		EXPECT_EQ(std::get<2>(minimizers.at(6)), 1);
		EXPECT_EQ(std::get<0>(minimizers.at(7)), 17);
		EXPECT_EQ(std::get<1>(minimizers.at(7)), 12);
		EXPECT_EQ(std::get<2>(minimizers.at(7)), 0);
		EXPECT_EQ(minimizers.size(), 8);
	}
}

// int main(int argc, char **argv) {
//   ::testing::InitGoogleTest(&argc, argv);
//   return RUN_ALL_TESTS();
// }