#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>
#include <functional>

#include "brown_alignment.hpp"

#define UP      'I'       // UP
#define LEFT    'D'       // LEFT
#define UPLEFT  'M'       // UP-LEFT
#define EMPTY   '0'       // EMPTY (LOCAL ALIGNMENT)

namespace brown {

typedef struct {
  int val;
  char parent;
} cell;

void print_val_matrix(cell** m, int rows, int cols) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      printf("%3d", m[i][j].val);
    }
    printf("\n");
  }
}

void print_char_matrix(cell** m, int rows, int cols) {
  printf("\n");
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      printf("  %c", m[i][j].parent);
    }
    printf("\n");
  }
}

// Matching function
int w(char q, char t, int match, int mismatch) {
  return (q == t) ? match : mismatch;
}

int maximum(int a, int b, int c, int d) {
  return std::max(std::max(a, b), std::max(c, d));
}

// Matrix initialization functions
void init_row(cell** m, int cols, int mult) {
  for (int i = 1; i < cols; ++i) {
    m[0][i].val = i * mult;
    m[0][i].parent = LEFT;
  }
}

void init_col(cell** m, int rows, int mult) {
  for (int i = 1; i < rows; ++i) {
    m[i][0].val = i * mult;
    m[i][0].parent = UP;
  }
}

void init_matrix(cell** m, int rows, int cols, int mult_row, int mult_col) {
  m[0][0].val = 0;
  m[0][0].parent = EMPTY;
  init_row(m, cols, mult_row);
  init_col(m, rows, mult_col);
}

// Alignment function for global alignment algorithms
void align_global(cell** m, int rows, int cols,
                  const char* query, const char* target,
                  int match, int mismatch, int gap) {

  for (int i = 1; i < rows; ++i) {
    for (int j = 1; j < cols; ++j) {
      int matched = m[i-1][j-1].val + w(query[i-1], target[j-1], match, mismatch);
      int insertion = m[i][j-1].val + gap;
      int deletion = m[i-1][j].val + gap;

      m[i][j].val = maximum(matched, insertion, deletion, std::numeric_limits<int>::min());

      if (m[i][j].val == matched) {
        m[i][j].parent = UPLEFT;
      } else if (m[i][j].val == insertion) {
        m[i][j].parent = LEFT;
      } else {
        m[i][j].parent = UP;
      }
    }
  }
}

// Alignment function for local alignment algorithms
// Stores the max value coordinates to prevent reiteration
int align_local(cell** m, int rows, int cols,
                const char* query, const char* target,
                int match, int mismatch, int gap,
                unsigned int& row, unsigned int& col) {

  int max = 0;

  for (int i = 1; i < rows; ++i) {
    for (int j = 1; j < cols; ++j) {
      int matched = m[i-1][j-1].val + w(query[i-1], target[j-1], match, mismatch);
      int insertion = m[i][j-1].val + gap;
      int deletion = m[i-1][j].val + gap;

      m[i][j].val = maximum(matched, insertion, deletion, 0);

      if (m[i][j].val > max) {
        max = m[i][j].val;
        row = i;
        col = j;
      }

      if (m[i][j].val == 0) {
        m[i][j].parent = EMPTY;
      } else if (m[i][j].val == matched && m[i-1][j-1].val != 0) {
        m[i][j].parent = UPLEFT;
      } else if (m[i][j].val == insertion && m[i][j-1].val != 0) {
        m[i][j].parent = LEFT;
      } else if(m[i][j].val == deletion && m[i-1][j].val != 0){
        m[i][j].parent = UP;
      } else{
        m[i][j].parent = EMPTY;
      }
    }
  }
  return max;
}

// Alignment score finder for last row and column
// Used for prefix-suffix, can be modified to fit other cases
int find_alignment_score(cell** m, int rows, int cols,
                         unsigned int& row, unsigned int& col) {

  row = rows-1;
  col = cols-1;
  int max = m[row][col].val;

  for (int i = rows-2; i >= 0; --i) {
    if (m[i][cols-1].val > max) {
      max = m[i][cols-1].val;
      row = i;
      col = cols-1;
    }
  }
  for (int i = cols-2; i >= 0; --i) {
    if (m[rows-1][i].val > max) {
      max = m[rows-1][i].val;
      row = rows-1;
      col = i;
    }
  }

  return m[row][col].val;
}

// CIGAR string finder
void find_cigar(cell** m, std::string& cigar, unsigned int& target_begin,
                int row, int col, std::function<bool(cell**, int&, int)> condition) {

  int counter = 0;
  char current = m[row][col].parent;
  target_begin = 0;

  while (condition(m, row, col)) {
    if (current != m[row][col].parent) {
      cigar.push_back(current);
      std::string num = std::to_string(counter);
      std::reverse(num.begin(), num.end());
      cigar.append(num);
      current = m[row][col].parent;
      counter = 1;
    } else {
      counter++;
    }
    if (m[row][col].parent == UPLEFT) {
      row--;
      col--;
    } else if (m[row][col].parent == UP) {
      row--;
    } else {
      col--;
    }
  }
  cigar.push_back(current);
  std::string num = std::to_string(counter);
  std::reverse(num.begin(), num.end());
  cigar.append(num);
  std::reverse(cigar.begin(), cigar.end());

  if(col != 0){
    target_begin = col;
  }

}

// Needleman-Wunsch algorithm for global sequence alignment
// Without CIGAR string
int needleman_wunsch(const char* query, int rows,
                     const char* target, int cols,
                     cell** m,
                     int match,
                     int mismatch,
                     int gap) {

  init_matrix(m, rows, cols, gap, gap);
  align_global(m, rows, cols, query, target, match, mismatch, gap);
  return m[rows-1][cols-1].val;
}

// With CIGAR string
int needleman_wunsch(const char* query, int rows,
                     const char* target, int cols,
                     cell** m,
                     int match,
                     int mismatch,
                     int gap,
                     std::string& cigar,
                     unsigned int& target_begin) {

  init_matrix(m, rows, cols, gap, gap);
  align_global(m, rows, cols, query, target, match, mismatch, gap);
  int alignment_score = m[rows-1][cols-1].val;

  print_val_matrix(m, rows, cols);
  print_char_matrix(m, rows, cols);

  find_cigar(m, cigar, target_begin, rows-1, cols-1, [](cell** m, int& row, int col)
    { return row > 0 || col > 0; });

  return alignment_score;
}

// Glocal (semi-global) alignment method for partial sequence alignment
// Without CIGAR string
int prefix_suffix(const char* query, int rows,
                  const char* target, int cols,
                  cell** m,
                  int match,
                  int mismatch,
                  int gap) {

  init_matrix(m, rows, cols, 0, 0);
  align_global(m, rows, cols, query, target, match, mismatch, gap);
  unsigned int row;
  unsigned int col;
  return find_alignment_score(m, rows, cols, row, col);
}

// With CIGAR string
int prefix_suffix(const char* query, int rows,
                  const char* target, int cols,
                  cell** m,
                  int match,
                  int mismatch,
                  int gap,
                  std::string& cigar,
                  unsigned int& target_begin) {

  init_matrix(m, rows, cols, 0, 0);
  align_global(m, rows, cols, query, target, match, mismatch, gap);
  unsigned int row;
  unsigned int col;
  int alignment_score = find_alignment_score(m, rows, cols, row, col);

  print_val_matrix(m, rows, cols);
  print_char_matrix(m, rows, cols);

  find_cigar(m, cigar, target_begin, row, col, [](cell** m, int& row, int col)
    { return row > 0; });
  return alignment_score;
}

// Smith-Waterman algorithm for local sequence alignment
// Without CIGAR string
int smith_waterman(const char* query, int rows,
                   const char* target, int cols,
                   cell** m,
                   int match,
                   int mismatch,
                   int gap) {

  init_matrix(m, rows, cols, 0, 0);
  unsigned int row;
  unsigned int col;
  align_local(m, rows, cols, query, target, match, mismatch, gap, row, col);
  return m[row][col].val;
}

// With CIGAR string
int smith_waterman(const char* query, int rows,
                   const char* target, int cols,
                   cell** m,
                   int match,
                   int mismatch,
                   int gap,
                   std::string& cigar,
                   unsigned int& target_begin) {

  init_matrix(m, rows, cols, 0, 0);
  unsigned int row;
  unsigned int col;
  align_local(m, rows, cols, query, target, match, mismatch, gap, row, col);
  int alignment_score = m[row][col].val;

  print_val_matrix(m, rows, cols);
  print_char_matrix(m, rows, cols);

  find_cigar(m, cigar, target_begin, row, col, [](cell** m, int& row, int col)
    { if (m[row][col].parent == EMPTY) {
        row = 0;
      }
      return row > 0 && col > 0;
    });
  return alignment_score;
}

// Driver function for multiple sequence alignment algorithms
// Without CIGAR string
int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match,
                       int mismatch,
                       int gap) {

  const int rows = query_length + 1;
  const int cols = target_length + 1;

  cell** m = new cell*[rows];
  for (int i = 0; i < rows; ++i) {
    m[i] = new cell[cols];
  }

  int alignment_score = 0;
  if (type == AlignmentType::global) {
    alignment_score = needleman_wunsch(query, rows, target, cols, m,
                                       match, mismatch, gap);
  } else if (type == AlignmentType::semi_global) {
    alignment_score = prefix_suffix(query, rows, target, cols, m,
                                    match, mismatch, gap);
  } else {
    alignment_score = smith_waterman(query, rows, target, cols, m,
                                     match, mismatch, gap);
  }

  for (int i = 0; i < rows; ++i) {
    delete m[i];
  }
  delete m;

  return alignment_score;
  }

  // With CIGAR string
int pairwise_alignment(const char* query, unsigned int query_length,
                       const char* target, unsigned int target_length,
                       AlignmentType type,
                       int match,
                       int mismatch,
                       int gap,
                       std::string& cigar,
                       unsigned int& target_begin) {

  const int rows = query_length + 1;
  const int cols = target_length + 1;

  cell** m = new cell*[rows];
  for (int i = 0; i < rows; ++i) {
    m[i] = new cell[cols];
  }

  int alignment_score = 0;
  if (type == AlignmentType::global) {
    alignment_score = needleman_wunsch(query, rows, target, cols, m,
                                       match, mismatch, gap, cigar, target_begin);
  } else if (type == AlignmentType::semi_global) {
    alignment_score = prefix_suffix(query, rows, target, cols, m,
                                    match, mismatch, gap, cigar, target_begin);
  } else {
    alignment_score = smith_waterman(query, rows, target, cols, m,
                                     match, mismatch, gap, cigar, target_begin);
  }

  for (int i = 0; i < rows; ++i) {
    delete m[i];
  }
  delete m;

  return alignment_score;
}

}
