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

	void middleMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length,  std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length, unsigned int border) { 

		long int position_of_last_found_min_mer = -1;
		unsigned int last_min_mer;

		std::unordered_map<unsigned int, std::pair<unsigned int, bool>> temp_results_cache;

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

				if(temp_results_cache.count(current_pos) != 0) {
					is_complement = temp_results_cache[current_pos].second;
					temp_mer = temp_results_cache[current_pos].first;
				} else {
					is_complement = false;
					temp_mer = convertKmerToInteger(temp_seg_ptr, j, k);

					if(temp_mer > border) {
						temp_mer = convertKmerToInteger(temp_seg_complement_ptr, j, k);
						is_complement = true;
					}

					temp_results_cache.emplace(current_pos, std::make_pair(temp_mer, is_complement));
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

			temp_results_cache.erase(i);
		}
	}

	void endFirstMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length, std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length, unsigned int border) {
		unsigned int min_mer;

		unsigned int stop_pos = std::min(std::get<1>(minimizers_vec[0]), window_length - 1);

		bool first_found = false;
		for(unsigned int i = 0; i < stop_pos; ++i) {
			bool is_complement = false;
			unsigned int temp_mer = convertKmerToInteger(sequence, i, k);

			if(temp_mer > border) {
				temp_mer = convertKmerToInteger(sequence_complement, i, k);
				is_complement = true;
			}

			if(!first_found || min_mer > temp_mer) {
				if(!first_found) first_found = true;

				min_mer = temp_mer;
				minimizers_vec.emplace_back(temp_mer, i, is_complement);
			}
		}
	}

	void endLastMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length, std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length, unsigned int border) {
		unsigned int min_mer;

		unsigned int stop_pos = std::max(std::get<1>(minimizers_vec[minimizers_vec.size() - 1]), sequence_length - k - window_length + 1);

		bool first_found = false;
		for(unsigned int i = sequence_length - k; i > stop_pos; --i) {
			bool is_complement = false;
			unsigned int temp_mer = convertKmerToInteger(sequence, i, k);

			if(temp_mer > border) {
				temp_mer = convertKmerToInteger(sequence_complement, i, k);
				is_complement = true;
			}

			if(!first_found || min_mer > temp_mer) {
				if(!first_found) first_found = true;

				min_mer = temp_mer;
				minimizers_vec.emplace_back(temp_mer, i, is_complement);
			}
		}
	}

	void endMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length, std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length, bool isStart, unsigned int border) {
		if(isStart) {
			endFirstMinimizers(sequence, sequence_complement, sequence_length, minimizers_vec, k, window_length, border);
		} else {	
			endLastMinimizers(sequence, sequence_complement, sequence_length, minimizers_vec, k, window_length, border);
		}		
	}

	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length) {
		
		std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vec;

		unsigned int border = (1 << (2 * k)) / 2;

		std::string result;		

		for(int i = 0; i < sequence_length; ++i) {
			result += complements_map[sequence[i]];
		}
	
		const char* seq_complement = result.c_str();

		middleMinimizers(sequence, seq_complement, sequence_length, minimizers_vec, k, window_length, border);
		
		endMinimizers(sequence, seq_complement, sequence_length, minimizers_vec, k, window_length, false, border);

		endMinimizers(sequence, seq_complement, sequence_length, minimizers_vec, k, window_length, true, border);

		return minimizers_vec;
	}
}
