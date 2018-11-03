#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>

#define MATCH 0
#define INSERT 1
#define DELETE 2

namespace orange {

    enum class AlignmentType {global, local, semi_global};

    namespace {
        
        typedef struct {
            int cost;
            int parent;
        } cell;

        void initMatrix (std::vector<std::vector<cell> > &matrix, int query_length, int target_length, int cost, AlignmentType type) {
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

            int i;
            for (i=1; i<=target_length; i++){
                matrix[0][i].cost = i*startRowMultiplyFactor;
		matrix[0][i].parent = INSERT;
            }
            for (i=1; i<=query_length; i++){
                matrix[i][0].cost = i*startColumnMultiplyFactor;
		matrix[i][0].parent = DELETE;
            }
        }

        int populateMatrix (std::vector<std::vector<cell> > &matrix, const char* query, const char* target, int match, int mismatch, int gap, AlignmentType type){
            int i, j, k;
            int costs[3];

	    int max = (type == AlignmentType::local) ? 0 : matrix[0][strlen(target)].cost;

            for (i=1; i <= strlen(query); i++) {
                for (j=1; j <= strlen(target); j++) {
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
			    max = std::max(max, matrix[i][j].cost);
			 }
		    } else if(type == AlignmentType::semi_global && j == strlen(target)) {
			max = std::max(max, matrix[i][j].cost);
		    }
                }
            }

	    return ((type == AlignmentType::local || type == AlignmentType::semi_global) ? max : matrix[strlen(query)][strlen(target)].cost);
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
			return populateMatrix(matrix, query, target, match , mismatch, gap, type);
                       }

    int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match,
                       int mismatch,
                       int gap,
                       std::string& cigar,
                       unsigned int& target_begin) {
                           
                           //TO DO
                           return 0;

                       }

}
