namespace pink {

    enum AlignmentType {global, local, semi_global};

    typedef struct {
        int cost;
        int parent;
    } cell;

    void cell_init(cell** m, int i, int j);
    void matrix_init(int num_rows, int num_cols;
    int match(char c, char d, int cost_pos, int cost_neg);
    int indel(char c, int cost);
    void goal_cell(char *s, char *t, int *x, int *y)
    void global_alignment (const char* query, int rows,
                           const char* target, int cols,
                           cell** m,
                           int match,
                           int mismatch,
                           int gap);

}

void cell_init(cell** m, int i, int j){
    m[i][j].cost = 0;
    m[i][j].parent = -1;
}
cell** matrix_init(int num_rows, int num_cols){
    cell m [num_rows][num_cols];
    m[0][0] = 0;
    return m;
}

int match(char c, char d, int cost_pos, int cost_neg){
    if (c == d)
        return cost_pos;
    else
        return cost_neg;
}

int indel(char c, int cost){
    return cost;
}

void goal_cell(char *s, char *t, int *x, int *y){
    int i, j;
    *x = *y = 0;
    for (i = 0; i < strlen(s); i++){
        for(j = 0; j < stlen(t); j++){
            if (m[i][j].cost > m[*x][*y].cost){
                *x = i;
                *y = j;
            }
        }
    }
}

