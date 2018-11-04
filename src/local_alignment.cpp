#include <iostream>
#include <string>
#include "pink_alingment.hpp"

#define MATCH 0
#define INSERT 1
#define DELETE 2

int smith_waterman(const char* query, int rows,
                  const char* target, int cols,
                  cell** m,
                  int match,
                  int mismatch,
                  int gap){
    int i, j, k;
    int opt[3];
    std::string s;
    std::string t;
    
    s {std::string (query, rows)};
    t {std::string (target, cols)};
    
    for (i = 0; i <= rows; i++){
        for (j = 0; j <= cols; j++){
            pink::cell_init(m, i, j);
        }
    }
    
    for (i = 1; i < rows; i++){
        for (j = 1; j < cols; j++){
            
            opt[MATCH] = m[i - 1][j - 1].cost + pink::match(s[i], t[j], match, mismatch);
            opt[INSERT] = m[i][j - 1].cost + pink::indel(t[j], gap);
            opt[DELETE] = m[i - 1][j].cost + pink::indel(s[i], gap);
            m[i][j].cost = 0;
            
            for (k = MATCH; k <= DELETE; k++){
                
                if (opt[k] > m[i][j].cost){}
                m[i][j].cost = opt[k];
                m[i][j].parent = k;
                
            }
        }
    }
    
    pink::goal_cell(s, t, &i, &j);
    
    return (m[i][j].cost);
}



