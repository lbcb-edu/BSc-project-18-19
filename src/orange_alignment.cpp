#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <list>
#include <utility>

#define STOP -1
#define MATCH 0
#define INSERT 1
#define DELETE 2

namespace orange {

	enum class AlignmentType {global, local, semi_global, no_alignment,};
        
        typedef struct {
            int cost;
            int parent;
        } cell;

	std::unordered_map<int, char> CIGAR_map = {
			{MATCH, 'M'},
			{INSERT, 'I'},
			{DELETE, 'D'}
	};

	void updatePosition(unsigned int& i, unsigned int& j, int parent) {
		switch(parent) {
			case MATCH:
				i = i - 1;
				j = j - 1;
				break;
			case INSERT:
				j = j - 1;
				break;
			case DELETE:
				i = i - 1;
				break;	
		}
	}

	std::string convertListToCIGARString(std::list<std::pair<char, unsigned int>> const &temp_list) {
		std::string CIGAR;

		for(auto const &pair : temp_list) {
			CIGAR += std::to_string(pair.second) + pair.first;
		}

		return CIGAR;
	}

	std::string constructCIGAR(std::vector<std::vector<cell>> const &matrix, unsigned int target_cell_row, unsigned int target_cell_column, unsigned int &target_begin, const char *query, const char *target, unsigned int query_length, unsigned int target_length) {
		unsigned int i = target_cell_row;
		unsigned int j = target_cell_column;

		bool firstIdentified = false;
		unsigned int counter;
		char lastChar, c;
		std::string prefix, suffix;

		if(i != query_length) 
			suffix = std::to_string(query_length - i) + 'S';
	
		std::list<std::pair<char, unsigned int>> temp_list;	

		while(matrix[i][j].parent != STOP) {
			if(matrix[i][j].parent==MATCH) {
				if(target[j-1]==query[i-1]) {
					c='=';
				}
				else {
					c='X';
				}
			} else 
				c = CIGAR_map.at(matrix[i][j].parent);

			if(firstIdentified && c == lastChar) {
				counter++;
			} else {
				if(firstIdentified) {
					temp_list.emplace_front(lastChar, counter);
				} else {
					firstIdentified = true;
				}

				lastChar = c;
				counter = 1;
			}

			updatePosition(i, j, matrix[i][j].parent);
		}

		if(i != 0)
			prefix = std::to_string(i) + 'S';

		target_begin=j;

		temp_list.emplace_front(lastChar, counter);

		return prefix + convertListToCIGARString(temp_list) + suffix;
	}

	void initMatrix (std::vector<std::vector<cell>> &matrix, int query_length, int target_length, int cost, AlignmentType type) {
		int startRowMultiplyFactor, startColumnMultiplyFactor;
		switch(type) {
			case AlignmentType::global :
				startRowMultiplyFactor = startColumnMultiplyFactor = cost;
				break;		
			case AlignmentType::local :
				startRowMultiplyFactor = startColumnMultiplyFactor = 0;
				break;
			case AlignmentType::semi_global :
				startRowMultiplyFactor = 0;
				startColumnMultiplyFactor = cost;
		}

		matrix[0][0].cost = 0;
		matrix[0][0].parent = STOP;

		int i;
		for (i=1; i<=target_length; ++i){
			matrix[0][i].cost = i*startRowMultiplyFactor;
			matrix[0][i].parent = (type == AlignmentType::global) ? INSERT : STOP;
            	}

            	for (i=1; i<=query_length; ++i){
                	matrix[i][0].cost = i*startColumnMultiplyFactor;
			matrix[i][0].parent = (type == AlignmentType::global) ? DELETE : STOP;
		}
	}

	int checkIfMaxAndUpdatePosition(int max, std::vector<std::vector<cell>> const &matrix, unsigned int i, unsigned int j, unsigned int& position_row, unsigned int& position_column) {
		int new_max = std::max(max, matrix[i][j].cost);
			
		if(new_max == matrix[i][j].cost) {
			position_row = i;
			position_column = j;
		}

		return new_max;
	}

	int populateMatrix (std::vector<std::vector<cell>> &matrix, const char* query, const char* target, unsigned int query_length, unsigned int target_length, int match, int mismatch, int gap, AlignmentType type, unsigned int& target_cell_row, unsigned int& target_cell_column){
		unsigned int i, j, k;
		int costs[3];

		int max = (type == AlignmentType::local) ? 0 : matrix[0][target_length].cost;

		target_cell_row = query_length;
		target_cell_column = target_length;

		for (i=1; i <= query_length; i++) {
			for (j=1; j <= target_length; j++) {
				costs[MATCH] = matrix[i-1][j-1].cost + (target[j-1]==query[i-1] ? match : mismatch);
				costs[INSERT] = matrix[i][j-1].cost + gap; 
				costs[DELETE] = matrix[i-1][j].cost + gap;

				matrix[i][j].cost = costs[MATCH];
				matrix[i][j].parent = MATCH;
				for (k=1; k<=2; k++) {
					if (costs[k] > matrix[i][j].cost) {
						matrix[i][j].cost = costs[k];
						matrix[i][j].parent = k;
                     			}
				}

				if(type == AlignmentType::local) {
					if(matrix[i][j].cost < 0) {
						matrix[i][j].cost = 0;
					} else {
						max = checkIfMaxAndUpdatePosition(max, matrix, i, j, target_cell_row, target_cell_column);
					} 

					if(matrix[i][j].cost == 0) {
						matrix[i][j].parent = STOP;
					}
				} else if(type == AlignmentType::semi_global && j == target_length) {
					max = checkIfMaxAndUpdatePosition(max, matrix, i, j, target_cell_row, target_cell_column);
				}
			}
		}

		switch(type) {
			case AlignmentType::local :
			case AlignmentType::semi_global :
					return max;
			default :
		   			return matrix[query_length][target_length].cost;
		}
	}

	int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match,
                       int mismatch,
                       int gap) {

		std::vector<std::vector<cell>> matrix (query_length+1, std::vector<cell>(target_length+1));

		initMatrix(matrix, query_length, target_length, gap, type);

		unsigned int temp_i;
		unsigned int temp_j;

		return populateMatrix(matrix, query, target, query_length, target_length, match , mismatch, gap, type, temp_i, temp_j);
	}

	int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match,
                       int mismatch,
                       int gap,
                       std::string& cigar,
                       unsigned int& target_begin) {
                           
		std::vector<std::vector<cell>> matrix (query_length+1, std::vector<cell>(target_length+1));

		initMatrix(matrix, query_length, target_length, gap, type);
			
		unsigned int target_cell_row;
		unsigned int target_cell_column;

		int score = populateMatrix(matrix, query, target, query_length, target_length, match , mismatch, gap, type, target_cell_row, target_cell_column);

		cigar = constructCIGAR(matrix, target_cell_row, target_cell_column, target_begin, query, target, query_length, target_length);

		return score;
	}

}
