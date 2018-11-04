#include <cstring>
#include <iostream>
//#include "pink_alignment.hpp"

namespace pink {
    
#define MATCH  0  
#define INSERT 1  
#define DELETE 2  
#define HEAD   -1  
    
typedef struct {
    int cost;  //cost of reaching this cell
    int parent;  //parent cell
} cell;

void matrix_init(cell** matrix, int rows, int columns, int gap) {
    for (int i = 0; i < rows; i++) {
        matrix[i][0].cost = i * gap;
        matrix[i][0].parent = HEAD;
    }
    
    for (int j = 1; j < columns; j++) {
        matrix[0][j].cost = j * gap;
        matrix[0][j] .parent = HEAD;
    }
}

int letter_match(char s, char t, int match, int mismatch) {
    if (s == t) {
        return match;
        
    } else {
        return mismatch;
    }
}

int string_compare(cell** matrix, const char *s, const char *t,
                                    unsigned rows, unsigned columns,
                                    int match, int mismatch, int gap) {
    unsigned i, j, k;  //counters
    int options[3];  //cost of options
    
    for (i = 1; i < rows; i++) {
        for (j = 1; j < columns; j++) {
            options[MATCH] = matrix[i-1][j-1].cost + letter_match(s[i - 1], t[j - 1], match, mismatch);
            options[INSERT] = matrix[i][j-1].cost + gap;
            options[DELETE] = matrix[i-1][j].cost + gap;
             
            matrix[i][j].cost = options[MATCH];
            matrix[i][j].parent = MATCH;
            
            for (k = 1; k <= 2; k++) {
                if (options[k] < matrix[i][j].cost) {
                    matrix[i][j].cost = options[k];
                    matrix[i][j].parent = k;
                }
            }
        }
    }
    
    std::cout << "\n~FINAL MATRIX~" << std::endl;
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < columns; j++) {
            std::cout << matrix[i][j].cost << " ";
        }
        std::cout << std::endl;
    }
    
    return matrix[rows - 1][columns - 1].cost;
}

int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       //AlignmentType type,
                       int match,
                       int mismatch,
                       int gap) {
    
    unsigned rows = query_length + 1;
    unsigned columns = target_length +1;
    
    cell** matrix = new cell*[rows];  //dynamical matrix
    for (unsigned i = 0; i < rows; i++) {
        matrix[i] = new cell[columns];
    }
    
    matrix_init(matrix, rows, columns, gap);
    
    std::cout << "\n~INITIAL MATRIX~" << std::endl;
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < columns; j++) {
            std::cout << matrix[i][j].cost << " ";
        }
        std::cout << std::endl;
    }
    
    return string_compare(matrix, query, target, rows, columns, match, mismatch, gap);
}
}

int main() {
    const char* query = "AAAC";
    const char* target = "AGC";
    
    int cost = pink::pairwise_alignment(query, strlen(query), target, strlen(target), 0, 1, 1);
                                
    std::cout << "\nFinal cost: " << cost << std::endl;
    return 0;
}

