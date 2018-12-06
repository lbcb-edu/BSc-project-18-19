#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>

#include "brown_minimizers.hpp"

namespace brown {

// C, A, T(U), G = 0, 1, 2, 3 for odd based
// C, A, T(U), G = 3, 2, 1, 0 for even based
std::map<char, int> odd_based;

// Position in array is even meaning that the actual position is odd
// e.g. array:  0 1 2 ...
//      actual: 1 2 3 ...
// Otherwise reverse odd-based
unsigned int value(const char* sequence, unsigned int pos, unsigned int k) {
  unsigned int val = 0;
  for (unsigned int i = 0; i < k; ++i) {
    if ((pos + i) % 2 == 0) {
      val = val * 10 + odd_based[sequence[pos + i]];
    } else {
      val = val * 10 + (3 - odd_based[sequence[pos + i]]);
    }
  }
  return val;
}

std::vector<std::pair<unsigned int, unsigned int>> minimizers(
    const char* sequence,
    unsigned int sequence_length,
    unsigned int k,
    unsigned int window_length) {

  odd_based['C'] = 0;
  odd_based['A'] = 1;
  odd_based['T'] = odd_based['U'] = 2;
  odd_based['G'] = 3;

  std::set<std::pair<unsigned int, unsigned int>> minimizers_set;
  std::pair<unsigned int, unsigned int> p;

  // Interior minimizers
  unsigned int l = window_length + k - 1;
  unsigned int last_window_pos = sequence_length - l;
  for (unsigned int i = 0; i <= last_window_pos; ++i) {
    unsigned int min_val = value(sequence, i, k);
    unsigned int min_pos = i;
    for (unsigned int j = 1; j < window_length; ++j) {
      unsigned int val = value(sequence, i + j, k);
      if (val < min_val) {
        min_val = val;
        min_pos = i + j;
      }
    }
    p = std::make_pair(min_val, min_pos);
    minimizers_set.insert(p);
  }

  // End minimizers
  // Front
  for (int u = 1; u < window_length; ++u) {
    unsigned int min_val = value(sequence, 0, k);
    unsigned int min_pos = 0;
    for (int i = 1; i < u; ++i) {
      unsigned int val = value(sequence, i, k);
      if (val < min_val) {
        min_val = val;
        min_pos = i;
      }
    }
    p = std::make_pair(min_val, min_pos);
    minimizers_set.insert(p);
  }
  // Back
  unsigned int last_pos = sequence_length - k;
  for (int u = 1; u < window_length; ++u) {
    unsigned int min_val = value(sequence, last_pos, k);
    unsigned int min_pos = last_pos;
    for (int i = 1; i < u; ++i) {
      unsigned int val = value(sequence, last_pos - i, k);
      if (val < min_val) {
        min_val = val;
        min_pos = last_pos - i;
      }
    }
    p = std::make_pair(min_val, min_pos);
    minimizers_set.insert(p);
  }

  std::vector<std::pair<unsigned int, unsigned int>> minimizers_vector(
      minimizers_set.begin(),
      minimizers_set.end());
  return minimizers_vector;
}

}
