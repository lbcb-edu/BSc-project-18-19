#ifndef BLUE_ALIGNMENT_H_INCLUDED
#define BLUE_ALIGNMENT_H_INCLUDED

#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <stdio.h>

enum AlignmentType {global, local, semi_global};

typedef struct {
    int value;
    std::pair<int, int> trace;
} Cell;

namespace blue
{
    void create_cigar_string(std::vector<std::vector<Cell> > &matrix, int i, int j, std::string &cigar, unsigned int &target_begin);

    void initialize_matrix(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           int match, int mismatch, int gap, std::vector<std::vector<Cell> > &matrix);

    int semi_global_alignment(const char* query, unsigned int query_length,
                                    const char* target, unsigned int target_length,
                                    int match, int mismatch, int gap, std::string &cigar, unsigned int &target_begin);

    int local_alignment(const char* query, unsigned int query_length,
                                 const char* target, unsigned int target_length,
                                 int match, int mismatch, int gap, std::string &cigar, unsigned int &target_begin);

    int global_alignment(const char* query, unsigned int query_length,
                                 const char* target, unsigned int target_length,
                                 int match, int mismatch, int gap, std::string &cigar, unsigned int &target_begin);

    int pairwise_alignment(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           AlignmentType type,
                           int match, int mismatch, int gap,
                           std::string &cigar, unsigned int &target_begin);

    int pairwise_alignment(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           AlignmentType type,
                           int match, int mismatch, int gap);

    AlignmentType getType(std::string type);
}

#endif // BLUE_ALIGNMENT_H_INCLUDED
