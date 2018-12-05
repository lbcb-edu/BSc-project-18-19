#include <iostream>
#include "pink_alignment.hpp"

namespace pink {
    
#define MATCH    0  
#define INSERT   1  
#define DELETE   2  
#define MISMATCH 3
#define HEAD    -1  
    
typedef struct {
    int cost;  
    int parent;  
} cell;

cell** matrix;
unsigned row;
unsigned col;

void init_row(int columns, int gap) {
    for (int j = 1; j < columns; j++) {
        matrix[0][j].cost = j * gap;
        matrix[0][j].parent = HEAD;
    }
}

void init_column(int rows, int gap) {
    for (int i = 0; i < rows; i++) {
        matrix[i][0].cost = i * gap;
        matrix[i][0].parent = HEAD;
    }
}

void init_matrix(AlignmentType type, int rows, int columns, int gap) {
    matrix[0][0].cost = 0;
    matrix[0][0].parent = HEAD;
    
    switch(type) {
        case global: init_row(columns, gap);
                     init_column(rows, gap);
                     break;
        
        //prefix-suffix
        case semi_global: init_row(columns, 0);
                          init_column(rows, gap);
                          break;
        
        case local: init_row(columns, 0);
                    init_column(rows, 0);
                    break;
    }
}

int compare_string(AlignmentType type,
                   const char *s, const char *t,
                   unsigned rows, unsigned columns,
                   int match, int mismatch, int gap) {
    unsigned i, j, k;  
    int options[3];  
    int max_cost = 0;
    row = 0;
    col = 0;
    
    if (type == semi_global) {
		row = 0;
        col = columns - 1;
    }
    
    for (i = 1; i < rows; i++) {
        for (j = 1; j < columns; j++) {
            int m = (s[i-1] == t[j-1]) ? match : mismatch;

            options[MATCH] = matrix[i-1][j-1].cost + m;
            options[INSERT] = matrix[i][j-1].cost + gap;
            options[DELETE] = matrix[i-1][j].cost + gap;
             
            matrix[i][j].cost = options[MATCH];
            matrix[i][j].parent = (m == match) ? MATCH : MISMATCH;
            
            for (k = 1; k <= 2; k++) {
                if (options[k] > matrix[i][j].cost) {
                    matrix[i][j].cost = options[k];
                    matrix[i][j].parent = k;
                }
            }
            
            switch (type) {
                case local: 
                    if (matrix[i][j].cost < 0) {
                        matrix[i][j].cost = 0;
                    } 
                    if (max_cost < matrix[i][j].cost) {   
                        max_cost = matrix[i][j].cost;
                        row = i;
                        col = j;
                    }
                    break;
                
                case semi_global:
                    if (j == columns - 1 && max_cost <= matrix[i][j].cost) {
                        max_cost = matrix[i][j].cost;
                        row = i;
                        col = j;
                    }
                    break;
                
                default: break;
            }
        }
    }
    
    if (type == global) {
        row = rows - 1;
        col = columns - 1;
        return matrix[rows - 1][columns - 1].cost;

    } else {
        return max_cost;
    }
}
    
char get_letter() { 
    switch (matrix[row][col].parent) {
        case MATCH:    return '=';
        case MISMATCH: return 'X';
        case INSERT:   return 'I';
        case DELETE:   return 'D';
        default:       return ' ';
    }
}
    
void switch_cell(char letter) {
    switch(letter) {
        case 'I': row--;
                  break;
        case 'D': col--;
                  break;
        default: row--;
                 col--;
                 break;
    }
} 
    
void make_cigar(AlignmentType type, 
			    std::string& cigar, unsigned int& target_begin) {
    
    switch(type) {
        case local: 
			while (matrix[row][col].cost != 0) {
            	char letter = get_letter();
            	cigar.push_back(letter);
            	switch_cell(letter);
        	}
			break;
        
        default: 
			while (matrix[row][col].parent != HEAD) {
	            char letter = get_letter();
	            cigar.push_back(letter);
	            switch_cell(letter);
	        }
        
	        while (row != 0) {
	            row--;
	            cigar.push_back('D');
	        }

	        while (col != 0) {
	            col--;
	            cigar.push_back('I');
	        }
	}
}

int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match, int mismatch, int gap) {

    unsigned rows = query_length + 1;
    unsigned columns = target_length +1;
    
    matrix = new cell*[rows];
    for (unsigned i = 0; i < rows; i++) {
        matrix[i] = new cell[columns];
    }
    
    init_matrix(type, rows, columns, gap);
    
    return compare_string(type, query, target, rows, columns, match, mismatch, gap);
}

int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match, int mismatch, int gap,
                       std::string& cigar, unsigned int& target_begin) {
    
    
    int cost = pairwise_alignment(query, query_length, target, target_length, type, match, mismatch, gap);
    
    make_cigar(type, cigar, target_begin);
    
    return cost;
}
}
