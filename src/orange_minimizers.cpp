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

	struct MerTupleHasher {
		std::size_t operator()(std::tuple<unsigned int, unsigned int, bool> const &val) const { return std::get<0>(val);}
	};

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

	void searchForMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length, std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec,
				unsigned int k, unsigned int window_length, unsigned int i_start, std::function<bool(unsigned int)> i_predicate,
				std::function<const char*(const char*, unsigned int)> segment_ptr_extractor, 
				std::function<int(unsigned int, bool)> j_start_calculator, std::function<bool(unsigned int, int)> j_predicate, std::function<int(int)> j_mover,
				std::function<unsigned int(unsigned int, int)> min_pos_calculator, unsigned int border, bool cache_cleaning_enabled) {
		
		long int position_of_last_found_min_mer = -1;
		unsigned int last_min_mer;

		std::unordered_map<unsigned int, std::pair<unsigned int, bool>> temp_results_cache;

		for(unsigned int i = i_start; i_predicate(i); ++i) {
			const char* temp_seg_ptr = segment_ptr_extractor(sequence, i);
			const char* temp_seg_complement_ptr = segment_ptr_extractor(sequence_complement, i);

			unsigned int min_mer = last_min_mer;
			bool is_min_complement;
			unsigned int min_pos = position_of_last_found_min_mer;
			bool first_found = position_of_last_found_min_mer >= i;

			for(int j = j_start_calculator(i, first_found); j_predicate(i, j); j = j_mover(j)) {
				unsigned int current_pos = min_pos_calculator(i, j);

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

			if(cache_cleaning_enabled) {
				temp_results_cache.erase(i);
			}
		}

	}

	void endMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length, std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length, bool isStart, unsigned int border) {
		unsigned int middle_min_pos;

		if(isStart) {
			middle_min_pos = std::get<1>(minimizers_vec[0]);
			searchForMinimizers(sequence, sequence_complement, sequence_length, minimizers_vec, k, window_length,
					k, [&](unsigned int i){ return i - k < middle_min_pos && i < k + window_length - 1;}, [](const char* seq, unsigned int i){ return seq;},
					[&](unsigned int i, bool condition) {return condition ? i - k : 0;}, [&](unsigned int i, int j){ return j + k <= i;},
					[&](int j){ return ++j;}, [](unsigned int i, int j){ return j;}, border, false);
		} else {
			middle_min_pos = std::get<1>(minimizers_vec[minimizers_vec.size() - 1]);	
			searchForMinimizers(sequence, sequence_complement, sequence_length, minimizers_vec, k, window_length,
					k, [&](unsigned int i){ return sequence_length - i > middle_min_pos && i < k + window_length - 1;}, 
					[&](const char* seq, unsigned int i){ return seq + sequence_length - i;},
					[&](unsigned int i, bool condition) {return condition ? 0 : i - k;}, [](unsigned int i, int j){ return j >= 0;},
					[&](int j){ return --j;}, [&](unsigned int i, int j){ return sequence_length - i + j;}, border, false);
		}		
	}

	void middleMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length,  std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length, unsigned int border) { 

		searchForMinimizers(sequence, sequence_complement, sequence_length, minimizers_vec, k, window_length,
					0, [&](unsigned int i){ return window_length + k + i - 1 <= sequence_length;}, [](const char* seq, unsigned int start){ return seq + start;},
					[&](unsigned int i, bool condition) {return condition ? window_length - 1 : 0;}, [&](unsigned int i, int j){ return j + k <= window_length + k - 1;},
					[](int j){ return ++j;}, [](unsigned int i, int j){ return i + j;}, border, true);
	}

	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length) {
		
		std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vec;

		unsigned int border = (1 << (2 * k)) / 2;

		const char* seq_complement = generateComplementSequence(sequence, sequence_length).c_str();

		middleMinimizers(sequence, seq_complement, sequence_length, minimizers_vec, k, window_length, border);

		endMinimizers(sequence, seq_complement, sequence_length, minimizers_vec, k, window_length, false, border);

		endMinimizers(sequence, seq_complement, sequence_length, minimizers_vec, k, window_length, true, border);

		return minimizers_vec;
	}
}
