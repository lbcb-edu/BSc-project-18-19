#include <map>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <tuple>
#include <algorithm>
#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <deque>

namespace orange {

	std::unordered_map<char, char> complements_map = {
		{'A', 'T'},
		{'T', 'A'},
		{'G', 'C'},
		{'C', 'G'}
	};

	std::unordered_map<char, int> values_map = {
		{'C', 0},
		{'A', 1},
		{'T', 2},
		{'G', 3}
	};

	std::string generateComplementSequence(const char *seq, unsigned int seq_length) {	
		std::string result;		

		for(int i = 0; i < seq_length; ++i) {
			result += complements_map[seq[i]];
		}

		return result;
	}

	unsigned int convertKmerToInteger(const char* temp_seg_ptr, int j, unsigned int k) {
		unsigned int result = 0;

		for(int i = 0; i < k; ++i) {
			unsigned int temp = (values_map[(temp_seg_ptr + j)[i]]) << (2 * (k - i - 1));
			result |= temp;
		}

		return result;
	}

	unsigned int convertKmerToIntegerWithPrevious(unsigned int last_kmer, char next_char, unsigned int mask) {
		return ((last_kmer << 2) | values_map[next_char]) & mask;
	}

	unsigned int convertKmerToIntegerWithPreviousInverse(unsigned int last_kmer, char next_char, unsigned int mask, unsigned int k) {
		return ((last_kmer >> 2) | (values_map[next_char] << 2*(k - 1))) & mask;
	}

	void middleMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length,  std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length, unsigned int mask) { 

		long int position_of_last_found_min_mer = -1;
		unsigned int last_min_mer;

		std::deque<std::tuple<unsigned int, unsigned int, unsigned int>> temp_results_cache;

		for(unsigned int i = 0; window_length + k + i - 1 <= sequence_length; ++i) {
			const char* temp_seg_ptr = sequence + i;
			const char* temp_seg_complement_ptr = sequence_complement + i;

			unsigned int min_mer = last_min_mer;
			bool is_min_complement;
			unsigned int min_pos = position_of_last_found_min_mer;
			bool first_found = position_of_last_found_min_mer >= i;

			for(int j = first_found ? window_length - 1 : 0; j + k <= window_length + k - 1; ++j) {
				unsigned int current_pos = i + j;;

				bool is_complement = false;
				unsigned int temp_mer;

				if(!temp_results_cache.empty() && current_pos <= std::get<2>(temp_results_cache[temp_results_cache.size() - 1])) {
					is_complement = std::get<1>(temp_results_cache[current_pos - i]) < std::get<0>(temp_results_cache[current_pos - i]);
					temp_mer = is_complement ? std::get<1>(temp_results_cache[current_pos - i]) : std::get<0>(temp_results_cache[current_pos - i]);
				} else {
					unsigned int comp;
					unsigned int not_comp;

					if(temp_results_cache.empty()) {
						not_comp = convertKmerToInteger(temp_seg_ptr, j, k);
						comp = convertKmerToInteger(temp_seg_complement_ptr, j, k);
					} else {
						not_comp = convertKmerToIntegerWithPrevious(std::get<0>(temp_results_cache[temp_results_cache.size() - 1]), *(sequence + current_pos + k - 1), mask);
						comp = convertKmerToIntegerWithPrevious(std::get<1>(temp_results_cache[temp_results_cache.size() - 1]), *(sequence_complement + current_pos + k - 1), mask);
					}

					is_complement = comp < not_comp;
					temp_mer = is_complement ? comp : not_comp;

					temp_results_cache.emplace_back(not_comp, comp, current_pos);
				}

				if(!first_found || min_mer > temp_mer) {
					if(!first_found) first_found = true;

					min_mer = temp_mer;
					is_min_complement = is_complement;
					min_pos = current_pos;
				}
			}
			
			if(min_pos != position_of_last_found_min_mer) {
				minimizers_vec.emplace_back(min_mer, min_pos, is_min_complement);
				position_of_last_found_min_mer = min_pos;
				last_min_mer = min_mer;
			}

			temp_results_cache.pop_front();
		}
	}

	void processEndMer(const char* sequence, const char* sequence_complement, std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec,unsigned int i, unsigned int k, unsigned int mask, unsigned int &min_mer, bool &first_found, unsigned int &last_not_comp, unsigned int &last_comp, bool isStart) {
		
		if(!first_found) {
			last_not_comp = convertKmerToInteger(sequence, i, k);
			last_comp = convertKmerToInteger(sequence_complement, i, k);
		} else {
			if(isStart) {
				last_not_comp = convertKmerToIntegerWithPrevious(last_not_comp, *(sequence + i + k - 1), mask);
				last_comp = convertKmerToIntegerWithPrevious(last_comp, *(sequence_complement + i + k - 1), mask);
			} else {
				last_not_comp = convertKmerToIntegerWithPreviousInverse(last_not_comp, *(sequence + i), mask, k);
				last_comp = convertKmerToIntegerWithPreviousInverse(last_comp, *(sequence_complement + i), mask, k);
			}
		}

		bool is_complement = last_comp < last_not_comp;
		unsigned int temp_mer = is_complement ? last_comp : last_not_comp;

		if(!first_found || min_mer > temp_mer) {
			if(!first_found) first_found = true;

			min_mer = temp_mer;
			minimizers_vec.emplace_back(temp_mer, i, is_complement);
		}
	}

	void endFirstMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length, std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length, unsigned int mask) {
		unsigned int min_mer;
		unsigned int stop_pos = std::min(std::get<1>(minimizers_vec[0]), window_length - 1);

		unsigned int last_not_comp;
		unsigned int last_comp;

		bool first_found = false;
		for(unsigned int i = 0; i < stop_pos; ++i) {
			processEndMer(sequence, sequence_complement, minimizers_vec,i, k, mask, min_mer, first_found, last_not_comp, last_comp, true);
		}
	}

	void endLastMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length, std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length, unsigned int mask) {
		unsigned int min_mer;

		unsigned int stop_pos = std::max(std::get<1>(minimizers_vec[minimizers_vec.size() - 1]), sequence_length - k - window_length + 1);

		unsigned int last_not_comp;
		unsigned int last_comp;

		bool first_found = false;
		for(unsigned int i = sequence_length - k; i > stop_pos; --i) {
			processEndMer(sequence, sequence_complement, minimizers_vec,i, k, mask, min_mer, first_found, last_not_comp, last_comp, false);
		}
	}

	void endMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length, std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length, bool isStart, unsigned int mask) {
		if(isStart) {
			endFirstMinimizers(sequence, sequence_complement, sequence_length, minimizers_vec, k, window_length, mask);
		} else {	
			endLastMinimizers(sequence, sequence_complement, sequence_length, minimizers_vec, k, window_length, mask);
		}		
	}

	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length) {
		
		std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vec;

		unsigned int mask = (1 << (2 * k)) - 1;

		std::string result;		

		for(int i = 0; i < sequence_length; ++i) {
			result += complements_map[sequence[i]];
		}
	
		const char* seq_complement = result.c_str();

		middleMinimizers(sequence, seq_complement, sequence_length, minimizers_vec, k, window_length, mask);
		
		endMinimizers(sequence, seq_complement, sequence_length, minimizers_vec, k, window_length, false, mask);

		endMinimizers(sequence, seq_complement, sequence_length, minimizers_vec, k, window_length, true, mask);

		return minimizers_vec;
	}
}
