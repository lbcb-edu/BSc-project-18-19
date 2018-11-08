#include "blue_alignment.hpp"
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <map>

namespace blue
{  
    typedef struct {
        int value;
        std::pair<int, int> trace;
    } Cell;
    
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
                           int match, int mismatch, int gap) {

                           std::string overloaded;
                           unsigned int broj;
                    return pairwise_alignment(query, query_length, target, target_length, type, match, mismatch, gap, overloaded, broj);
    }

    int pairwise_alignment(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           AlignmentType type,
                           int match, int mismatch, int gap,
                           std::string &cigar, unsigned int &target_begin) {

        switch (type) {
            case global:
                return global_alignment(query, query_length, target, target_length, match, mismatch, gap, cigar, target_begin);
            case local:
               return local_alignment(query, query_length, target, target_length, match, mismatch, gap, cigar, target_begin);

            case semi_global:
                return semi_global_alignment(query, query_length, target, target_length, match, mismatch, gap, cigar, target_begin);

            default:
                std::cout << "Wrong Alignment type!" << std::endl;
                break;
        }
        return -1;
    }
    //************************************************************************************************************************
    int local_alignment(const char* query, unsigned int query_length,
                                 const char* target, unsigned int target_length,
                                 int match, int mismatch, int gap,
                                 std::string &cigar, unsigned int &target_begin) {

        std::vector<std::vector<Cell> > matrix;
        matrix.resize(query_length+1, std::vector<Cell>(target_length+1));

        matrix[0][0].value = 0;
        matrix[0][0].trace.first = -1;
        matrix[0][0].trace.second = -1;

        for (int i=1; i < query_length + 1; i++){
            matrix[i][0].value = 0;
            matrix[i][0].trace.first = -1;
            matrix[i][0].trace.second = -1;
        }

        for (int i=1; i < target_length + 1; i++){
            matrix[0][i].value = 0;
            matrix[0][i].trace.first = -1;
            matrix[0][i].trace.second = -1;
        }

        std::pair<int, int> max_cell;
        max_cell.first = 0;
        max_cell.second = 0;
        int max_val = 0;

        for(int i = 1; i < query_length+1; i++) {
            for(int j = 1; j<target_length+1; j++ ) {
                int match_mismatch = query[i-1] == target[j-1] ?
                                        matrix[i-1][j-1].value + match : matrix[i-1][j-1].value + mismatch;
                int insertion = matrix[i][j-1].value + gap;
                int deletion = matrix[i-1][j].value + gap;

                if ( match_mismatch < 0 && insertion < 0 && deletion < 0) {
                    matrix[i][j].value = 0;
                    matrix[i][j].trace.first = -1;
                    matrix[i][j].trace.second = -1;
                } else {
                    if(insertion > deletion) {
                        if(match_mismatch > insertion) {
                            matrix[i][j].value = match_mismatch;
                            matrix[i][j].trace.first = i-1;
                            matrix[i][j].trace.second = j-1;

                        } else if (match_mismatch == insertion) {
                            if (matrix[i-1][j-1].value > matrix[i][j-1].value) {
                                matrix[i][j].trace.first = i-1;
                                matrix[i][j].trace.second = j-1;
                            } else {
                                matrix[i][j].trace.first = i;
                                matrix[i][j].trace.second = j-1;
                            }
                            matrix[i][j].value = match_mismatch;

                        } else {
                            matrix[i][j].value = insertion;
                            matrix[i][j].trace.first = i;
                            matrix[i][j].trace.second = j-1;
                        }

    //***********************************************************************************
                    } else {
                        if(match_mismatch > deletion) {
                            matrix[i][j].value = match_mismatch;
                            matrix[i][j].trace.first = i-1;
                            matrix[i][j].trace.second = j-1;

                        } else if (match_mismatch == deletion) {
                            if (matrix[i-1][j-1].value > matrix[i-1][j].value) {
                                matrix[i][j].trace.first = i-1;
                                matrix[i][j].trace.second = j-1;
                            } else {
                                matrix[i][j].trace.first = i-1;
                                matrix[i][j].trace.second = j;
                            }
                            matrix[i][j].value = match_mismatch;
                        } else {
                            matrix[i][j].value = deletion;
                            matrix[i][j].trace.first = i-1;
                            matrix[i][j].trace.second = j;
                        }
                    }
                    if (matrix[i][j].value > max_val) {
                        max_val = matrix[i][j].value;
                        max_cell.first = i;
                        max_cell.second = j;
                    }
                }
            }
        }

        create_cigar_string(matrix, max_cell.first, max_cell.second, cigar, target_begin);
        return max_val;
    }

    int global_alignment(const char* query, unsigned int query_length,
                                 const char* target, unsigned int target_length,
                                 int match, int mismatch, int gap,
                                 std::string &cigar, unsigned int &target_begin) {
        std::vector<std::vector<Cell> > matrix;
        matrix.resize(query_length+1, std::vector<Cell>(target_length+1));

        matrix[0][0].value = 0;
        matrix[0][0].trace.first = -1;
        matrix[0][0].trace.second = -1;

        for(int i = 1; i < query_length+1; i++) {
            matrix[i][0].value = gap*i;
            matrix[i][0].trace.first = i-1;
            matrix[i][0].trace.second = 0;
        }

        for(int i = 1; i<target_length + 1; i++) {
            matrix[0][i].value = gap*i;
            matrix[0][i].trace.first = 0;
            matrix[0][i].trace.second = i-1;
        }

        initialize_matrix(query, query_length, target, target_length, match, mismatch, gap, matrix);
        create_cigar_string(matrix, query_length, target_length, cigar, target_begin);
        return matrix[query_length][target_length].value;
    }

    int semi_global_alignment(const char* query, unsigned int query_length,
                                    const char* target, unsigned int target_length,
                                    int match, int mismatch, int gap,
                                    std::string &cigar, unsigned int &target_begin) {

        std::vector<std::vector<Cell> > matrix;
        matrix.resize(query_length+1, std::vector<Cell>(target_length+1));

        matrix[0][0].value = 0;
        matrix[0][0].trace.first = -1;
        matrix[0][0].trace.second = -1;

        for (int i=1; i < query_length + 1; i++){
            matrix[i][0].value = 0;
            matrix[i][0].trace.first = -1;
            matrix[i][0].trace.second = -1;
        }

        for (int i=1; i < target_length + 1; i++){
            matrix[0][i].value = 0;
            matrix[0][i].trace.first = -1;
            matrix[0][i].trace.second = -1;
        }

        initialize_matrix(query, query_length, target, target_length, match, mismatch, gap, matrix);

        std::pair<int, int> max_cell;
        max_cell.first = 0;
        max_cell.second = 0;
        int max_val = matrix[0][0].value;

        for (int i = 1; i < query_length + 1; i++) {
            if (matrix[i][target_length].value > max_val) {
                max_val = matrix[i][target_length].value;
                max_cell.first = i;
                max_cell.second = target_length;
            }
        }

        for (int i = 1; i < target_length + 1; i++) {
            if (matrix[query_length][i].value > max_val) {
                max_val = matrix[query_length][i].value;
                max_cell.first = query_length;
                max_cell.second = i;
            }
        }

        create_cigar_string(matrix, max_cell.first, max_cell.second, cigar, target_begin);
        return max_val;
    }
    //************************************************************************************************************************

    void initialize_matrix(const char* query, unsigned int query_length,
                           const char* target, unsigned int target_length,
                           int match, int mismatch, int gap, std::vector<std::vector<Cell> > &matrix) {

        for(int i = 1; i < query_length+1; i++) {
            for(int j = 1; j<target_length+1; j++ ) {
                int match_mismatch = query[i-1] == target[j-1] ?
                                        matrix[i-1][j-1].value + match : matrix[i-1][j-1].value + mismatch;
                int insertion = matrix[i][j-1].value + gap;
                int deletion = matrix[i-1][j].value + gap;


                if(insertion > deletion) {
                    if(match_mismatch > insertion) {
                        matrix[i][j].value = match_mismatch;
                        matrix[i][j].trace.first = i-1;
                        matrix[i][j].trace.second = j-1;

                    } else if (match_mismatch == insertion) {
                        if (matrix[i-1][j-1].value > matrix[i][j-1].value) {
                            matrix[i][j].trace.first = i-1;
                            matrix[i][j].trace.second = j-1;
                        } else {
                            matrix[i][j].trace.first = i;
                            matrix[i][j].trace.second = j-1;
                        }
                        matrix[i][j].value = match_mismatch;

                    } else {
                        matrix[i][j].value = insertion;
                        matrix[i][j].trace.first = i;
                        matrix[i][j].trace.second = j-1;
                    }

    //***********************************************************************************
                } else {  //deletion > insertion
                    if(match_mismatch > deletion) {
                        matrix[i][j].value = match_mismatch;
                        matrix[i][j].trace.first = i-1;
                        matrix[i][j].trace.second = j-1;

                    } else if (match_mismatch == deletion) {
                        if (matrix[i-1][j-1].value > matrix[i-1][j].value) {
                            matrix[i][j].trace.first = i-1;
                            matrix[i][j].trace.second = j-1;
                        } else {
                            matrix[i][j].trace.first = i-1;
                            matrix[i][j].trace.second = j;
                        }
                        matrix[i][j].value = match_mismatch;

                    } else {
                        matrix[i][j].value = deletion;
                        matrix[i][j].trace.first = i-1;
                        matrix[i][j].trace.second = j;
                    }
                }
            }
        }
    }

    void create_cigar_string(std::vector<std::vector<Cell> > &matrix, int start_row, int start_column,
                                    std::string &adress, unsigned int &target_begin) {
        int mis_counter, del_counter, ins_counter;
        mis_counter = del_counter = ins_counter = 0;

        std::string cigar_reverse;
        Cell current = matrix[start_row][start_column];

        std::pair<int, int> pos;
        pos.first = start_row;
        pos.second = start_column;

        while( current.trace.first != -1 && current.trace.second != -1) {
            //provjera za match
            if(current.trace.first == pos.first-1 && current.trace.second == pos.second-1) {
                if(del_counter != 0) {
                    cigar_reverse+=std::to_string(del_counter);
                    del_counter = 0;
                } else if (ins_counter != 0) {
                    cigar_reverse +=std::to_string(ins_counter);
                    ins_counter = 0;
                }
                if(mis_counter == 0) {
                    cigar_reverse+='M';
                }
                ++mis_counter;

                // provjera za insertion
            } else if(current.trace.first == pos.first && current.trace.second == pos.second-1) {
                if(del_counter != 0) {
                    cigar_reverse+=std::to_string(del_counter);
                    del_counter = 0;
                } else if (mis_counter != 0) {
                    cigar_reverse +=std::to_string(mis_counter);
                    mis_counter = 0;
                }
                if(ins_counter == 0) {
                    cigar_reverse+='I';
                }
                ++ins_counter;

                //deletion
            } else {
                if(ins_counter != 0) {
                    cigar_reverse+=std::to_string(ins_counter);
                    ins_counter = 0;
                } else if (mis_counter != 0) {
                    cigar_reverse +=std::to_string(mis_counter);
                    mis_counter = 0;
                }
                if(del_counter == 0) {
                    cigar_reverse+='D';
                }
                ++del_counter;
            }
            pos = current.trace;
            current = matrix[current.trace.first][current.trace.second];
        }

       target_begin = (unsigned int) pos.second;

        if(mis_counter != 0) {
            cigar_reverse+=std::to_string(mis_counter);
        } else if(ins_counter != 0){
            cigar_reverse+=std::to_string(ins_counter);
        } else {
            cigar_reverse+=std::to_string(del_counter);
        }

        std::string number = "";
        std::string cigar = "";
        for(int i=cigar_reverse.length()-1; i>=0; i-- ) {
            char c = cigar_reverse.at(i);
            if(isdigit(c)) {
                number.insert(0, 1, c);
            } else {
                cigar += number;
                cigar += c;
                number = "";
            }
        }

        adress = cigar;
        return;
    }

    AlignmentType getType(std::string type) {
        if(type.compare("global") == 0) return global;
        if(type.compare("local") == 0) return local;
        if(type.compare("semi_global")==0) return semi_global;
    }
}
