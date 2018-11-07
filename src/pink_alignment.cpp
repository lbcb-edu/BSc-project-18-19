#include <iostream>
#include "pink_alignment.hpp"

namespace pink {
    
#define MATCH   0  
#define INSERT  1  
#define DELETE  2  
#define HEAD   -1  
    
typedef struct {
    int cost;  
    int parent;  
} cell;

void row_init(cell** matrix, int columns, int gap) {
    for (int j = 1; j < columns; j++) {
        matrix[0][j].cost = j * gap;
        matrix[0][j].parent = HEAD;
    }
}

void column_init(cell** matrix, int rows, int gap) {
    for (int i = 0; i < rows; i++) {
        matrix[i][0].cost = i * gap;
        matrix[i][0].parent = HEAD;
    }
}

void matrix_init(cell** matrix, AlignmentType type, int rows, int columns, int gap) {
    matrix[0][0].cost = 0;
    matrix[0][0].parent = HEAD;
    
    switch(type) {
        case global: row_init(matrix, columns, gap);
                     column_init(matrix, rows, gap);
                     break;
        
        //prefix-suffix
        case semi_global: row_init(matrix, columns, 0);
                          column_init(matrix, rows, gap);
                          break;
        
        case local: row_init(matrix, columns, 0);
                    column_init(matrix, rows, 0);
                    break;
    }
}

int letter_match(char s, char t, int match, int mismatch) {
    if (s == t) {
        return match;
        
    } else {
        return mismatch;
    }
}

int string_compare(cell** matrix, AlignmentType type,
                   const char *s, const char *t,
                   unsigned rows, unsigned columns,
                   int match, int mismatch, int gap,
                   unsigned* row, unsigned* col) {
    unsigned i, j, k;  
    int options[3];  
    int max_cost = 0;
    
    if (type == semi_global) {
		*row = 0;
        *col = columns - 1;
    }
    
    for (i = 1; i < rows; i++) {
        for (j = 1; j < columns; j++) {
            options[MATCH] = matrix[i-1][j-1].cost + letter_match(s[i - 1], t[j - 1], match, mismatch);
            options[INSERT] = matrix[i][j-1].cost + gap;
            options[DELETE] = matrix[i-1][j].cost + gap;
             
            matrix[i][j].cost = options[MATCH];
            matrix[i][j].parent = MATCH;
            
            for (k = 1; k <= 2; k++) {
                if (options[k] > matrix[i][j].cost) {
                    matrix[i][j].cost = options[k];
                    matrix[i][j].parent = k;
                }
            }
            
            if (type == local) {
                if (matrix[i][j].cost < 0) {
                    matrix[i][j].cost = 0;
                } 
                if (max_cost < matrix[i][j].cost) {   
                    max_cost = matrix[i][j].cost;
                    *row = i;
                    *col = j;
                }
            }
            
            if (type == semi_global && j == columns - 1 && max_cost <= matrix[i][j].cost) {
                max_cost = matrix[i][j].cost;
                *row = i;
                *col = j;
            }
        }
    }
    
//    std::cout << "\n~FINAL MATRIX~" << std::endl;
//    for (unsigned i = 0; i < rows; i++) {
//        for (unsigned j = 0; j < columns; j++) {
//            std::cout << matrix[i][j].cost << " ";
//        }
//        std::cout << std::endl;
//    }
    
    if (type == global) {
        return matrix[rows - 1][columns - 1].cost;

    } else {
        return max_cost;
    }
}
    
char switch_cell (unsigned* i, unsigned* j, cell** matrix) {
	char position;
        
    switch (matrix[*i][*j].parent) {
        case MATCH:  (*i)--;
               		 (*j)--;
               		 position = 'M';
               		 break;

        case INSERT: (*j)--;
               		 position = 'I';
               		 break;

        case DELETE: (*i)--;
               		 position = 'D';
               		 break;
    }
    
	return position;
}
    
    
void makeCigar(cell** matrix, AlignmentType type, unsigned row, unsigned col, 
			   std::string& cigar, unsigned int& target_begin) {
    unsigned i = row;
    unsigned j = col;
   
    if (type == global || type == semi_global) {
        while (matrix[i][j].parent != HEAD) {
            cigar.push_back(switch_cell(&i, &j, matrix));
        }
        
        while (i != 0) {
            i--;
            cigar.push_back('D');
        }

        while (j != 0) {
            j--;
            cigar.push_back('I');
        }

    } else if (type == local) {
        while (matrix[i][j].cost != 0){
             cigar.push_back(switch_cell(&i, &j, matrix));
        }
    }
        
    if (col != 0) {
        target_begin = col;
    }
}

int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match, int mismatch, int gap) {

    unsigned rows = query_length + 1;
    unsigned columns = target_length +1;
    
    cell** matrix = new cell*[rows];
    for (unsigned i = 0; i < rows; i++) {
        matrix[i] = new cell[columns];
    }
    
    matrix_init(matrix, type, rows, columns, gap);
    
//    std::cout << "\n~INITIAL MATRIX~" << std::endl;
//    for (unsigned i = 0; i < rows; i++) {
//        for (unsigned j = 0; j < columns; j++) {
//            std::cout << matrix[i][j].cost << " ";
//        }
//        std::cout << std::endl;
//    }
    
    unsigned row;
    unsigned col;
    
    return string_compare(matrix, type, query, target, rows, columns, match, mismatch, gap, &row, &col);
}
    

int pairwise_alignment(const char* query, unsigned int query_length,
                        const char* target, unsigned int target_length,
                        AlignmentType type,
                        int match, int mismatch, int gap,
                        std::string& cigar, unsigned int& target_begin){
    
    unsigned rows = query_length + 1;
    unsigned columns = target_length + 1;
    
    cell** matrix = new cell*[rows];
    for (unsigned i = 0; i < rows; i++){
        matrix[i] = new cell[columns];
    }
    
    matrix_init(matrix, type, rows, columns, gap);
    
    unsigned row;
    unsigned col;
    
    int cost = string_compare(matrix, type, query, target, rows, columns, match, mismatch, gap, &row, &col);
    
    makeCigar(matrix, type, row, col, cigar, target_begin);
    
    return cost;
}
}
