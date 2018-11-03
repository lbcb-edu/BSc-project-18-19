#include <iostream>
#include <string>
#define MATCH 0
#define INSERT 1
#define DELETE 2
#include "pink_alingment.hpp"


int needleman_wunsch(const char* query, int rows,
                     const char* target, int cols,
                     cell** m,
                     int match,
                     int mismatch,
                     int gap){
    
    int i, j;
    int opt[3];
    int d = 1;
    std::string s;
    std::string t;
    
    s = std::string(query, rows);
    t = std::string(target, cols);
    
    int max;
    if (rows > cols){
        max = rows;
    } else {
        max = cols;
    }
    
    for (i = 0; i < max; i++) {
        m[0][i].cost = i * d;
        m[i][0].cost = i * d;
    }
    
    for (i = 1; i < rows; ++i) {
        for (j = 1; j < cols; ++j) {
            
            opt[MATCH] = m[i - 1][j - 1].cost + pink::match(query[i-1], target[j-1], match, mismatch);
            opt[INSERT] = m[i][j - 1].cost + pink::indel(t[j], gap);
            opt[DELETE] = m[i - 1][j].cost + pink::indel(s[i], gap);
           
            m[i][j].cost = 0;
            
            for (k = MATCH; k <= DELETE; k++){
                
                if (opt[k] > m[i][j].cost){
                    m[i][j].cost = opt[k];
                    m[i][j].parent = k;
                
                }
            }
        }
    }
    
    pink::goal_cell(rows, cols, &i, &j);
    
    return m[i][j].cost;
}
