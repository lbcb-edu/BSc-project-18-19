#include <stdlib.h>
#include <string>

namespace orange {

    enum class AlignmentType {global, local, semi_global};

    int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match,
                       int mismatch,
                       int gap) {

                           //TO DO
                           return 0;

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