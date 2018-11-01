#include <stdlib.h>
#include <string>
#include <vector>

#define MATCH 0
#define INSERT 1
#define DELETE 2

namespace orange {

    namespace {
        
        typedef struct {
            int cost;
            int parent;
        } cell;

        void initMatrix (std::vector<std::vector<cell> > &matrix, int query_length, int target_length, int cost) {
            int i;
            for (i=0; i<=target_length; i++){
                matrix[0][i].cost = i*cost;
            }
            for (i=0; i<=query_length; i++){
                matrix[i][0].cost = i*cost;
            }
        }

        void populateMatrix (std::vector<std::vector<cell> > &matrix, const char* query, const char* target, int match, int mismatch, int gap){
            int i, j, k;
            int costs[3];
            for (i=1; i <= strlen(query); i++) {
                for (j=1; j <= strlen(target); j++) {
                    costs[MATCH] = matrix[i-1][j-1].cost + (target[j-1]==query[i-1] ? match : mismatch);
                    costs[INSERT] = matrix[i][j-1].cost + gap; 
                    costs[DELETE] = matrix[i-1][j].cost + gap;
                    matrix[i][j].cost = costs[MATCH];
                    matrix[i][j].parent = MATCH;
                    for (k=1; k<=2; k++) {
                        if (costs[k] < matrix[i][j].cost) {
                            matrix[i][j].cost = costs[k];
                            matrix[i][j].parent = k;
                        }
                    }
                }
            }
        }


    }

    enum class AlignmentType {global, local, semi_global};

    int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match,
                       int mismatch,
                       int gap) {
                        std::vector<std::vector<cell>> matrix (query_length+1, std::vector<cell>(target_length+1));
                        switch(type) {
                            case AlignmentType::global:
                                initMatrix(matrix, query_length, target_length, gap);
                                populateMatrix(matrix, query, target, match , mismatch, gap);
                                return (matrix[query_length][target_length].cost);
                            case AlignmentType::semi_global:
                                
                                //TO DO

                                return 0;
                            case AlignmentType::local:

                                //TO DO

                                return 0;
                        }

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