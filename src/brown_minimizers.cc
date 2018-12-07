#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <tuple>

#include "brown_minimizers.hpp"

namespace brown {

// C, A, T/U, G = 0, 1, 2, 3 for odd based
// C, A, T/U, G = 3, 2, 1, 0 for even based
// Position in array is even meaning that the actual position is odd
// e.g. array:  0 1 2 ...
//      actual: 1 2 3 ...
// Otherwise reverse odd-based
std::map<char, int> char_to_val = {{'C', 0}, {'A', 1}, {'T', 2}, {'U', 2}, {'G', 3}};

unsigned int value(const char* sequence, unsigned int pos, unsigned int k) {
  unsigned int val = 0;
  for (unsigned int i = 0; i < k; ++i) {
    val *= 10;
    if ((pos + i) % 2 == 0) {
      val += char_to_val[sequence[pos + i]];
    } else {
      val += (3 - char_to_val[sequence[pos + i]]);
    }
  }
  return val;
}

unsigned int value_reverse_complement(const char* sequence,
                                      unsigned int sequence_length,
                                      unsigned int pos,
                                      unsigned int k) {

  unsigned int val = 0;
  for (unsigned int i = 0; i < k; ++i) {
    val *= 10;
    if ((sequence_length - 1 - (pos - i)) % 2 == 0) {
      val += (3 - char_to_val[sequence[pos - i]]);
    } else {
      val += char_to_val[sequence[pos - i]];
    }
  }
  return val;
}

void interior_minimizers_fill(
    std::set<std::tuple<unsigned int, unsigned int, bool>>& minimizers_set,
    const char* sequence,
    unsigned int sequence_length,
    unsigned int k,
    unsigned int window_length) {

  std::tuple<unsigned int, unsigned int, bool> t;

  unsigned int l = window_length + k - 1;
  unsigned int last_window_pos = sequence_length - l;

  for (unsigned int i = 0; i <= last_window_pos; ++i) {
    unsigned int min_val = value(sequence, i, k);
    unsigned int min_pos = i;
    bool         min_strand = 0;
    for (unsigned int j = 0; j < window_length; ++j) {
      unsigned int val = value(sequence, i + j, k);
      unsigned int val_rev_com = value_reverse_complement(sequence, sequence_length, i + l - 1 - j, k);
      if (val < min_val) {
        min_val = val;
        min_pos = i + j;
        min_strand = 0;
      }
      if (val_rev_com < min_val) {
        min_val = val_rev_com;
        min_pos = sequence_length - 1 - (i + l - 1 - j);
        min_strand = 1;
      }
    }
    t = std::make_tuple(min_val, min_pos, min_strand);
    minimizers_set.insert(t);
  }
}

void front_end_minimizers_fill(
    std::set<std::tuple<unsigned int, unsigned int, bool>>& minimizers_set,
    const char* sequence,
    unsigned int sequence_length,
    unsigned int k,
    unsigned int window_length) {

  std::tuple<unsigned int, unsigned int, bool> t;

  for (int u = 1; u < window_length; ++u) {
    unsigned int min_val = value(sequence, 0, k);
    unsigned int min_pos = 0;
    bool         min_strand = 0;
    for (int i = 0; i < u; ++i) {
      unsigned int val = value(sequence, i, k);
      unsigned int val_rev_com = value_reverse_complement(sequence, sequence_length, u + k - 2 - i, k);
      if (val < min_val) {
        min_val = val;
        min_pos = i;
        min_strand = 0;
      }
      if (val_rev_com < min_val) {
        min_val = val_rev_com;
        min_pos = sequence_length - 1 - (u + k - 2 - i);
        min_strand = 1;
      }
    }
    t = std::make_tuple(min_val, min_pos, min_strand);
    minimizers_set.insert(t);
  }
}

void back_end_minimizers_fill(
    std::set<std::tuple<unsigned int, unsigned int, bool>>& minimizers_set,
    const char* sequence,
    unsigned int sequence_length,
    unsigned int k,
    unsigned int window_length) {

  std::tuple<unsigned int, unsigned int, bool> t;

  unsigned int last_pos = sequence_length - k;
  for (int u = 1; u < window_length; ++u) {
    unsigned int min_val = value(sequence, last_pos, k);
    unsigned int min_pos = last_pos;
    bool         min_strand = 0;
    for (int i = 1; i < u; ++i) {
      unsigned int val = value(sequence, last_pos - i, k);
      unsigned int val_rev_com = value_reverse_complement(sequence, sequence_length, sequence_length - 1 - i, k);
      if (val < min_val) {
        min_val = val;
        min_pos = last_pos - i;
        min_strand = 0;
      }
      if (val_rev_com < min_val) {
        min_val = val_rev_com;
        min_pos = sequence_length - 1 - i;
        min_strand = 1;
      }
    }
    t = std::make_tuple(min_val, min_pos, min_strand);
    minimizers_set.insert(t);
  }
}

// bool | 0 == original strand minimizer
//      | 1 == reverse complement strand minimizer
std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(
    const char* sequence,
    unsigned int sequence_length,
    unsigned int k,
    unsigned int window_length) {

  if (sequence_length < window_length + k - 1) {
    std::vector<std::tuple<unsigned int, unsigned int, bool>> empty_vector;
    return empty_vector;
  }

  std::set<std::tuple<unsigned int, unsigned int, bool>> minimizers_set;
  
  interior_minimizers_fill(minimizers_set, sequence, sequence_length, k, window_length);

  front_end_minimizers_fill(minimizers_set, sequence, sequence_length, k, window_length);

  back_end_minimizers_fill(minimizers_set, sequence, sequence_length, k, window_length);

  std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector(
      minimizers_set.begin(), minimizers_set.end());
  return minimizers_vector;
}

}