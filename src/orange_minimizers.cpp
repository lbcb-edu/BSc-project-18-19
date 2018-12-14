#include <map>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <tuple>
#include <algorithm>
#include <functional>
#include <unordered_set>

namespace orange {

	struct MerTupleHasher {
		std::size_t operator()(std::tuple<unsigned int, unsigned int, bool> const &val) const { return std::get<1>(val);}
	};

	std::map<char, char> complements_map = {
		{'A', 'T'},
		{'T', 'A'},
		{'G', 'C'},
		{'C', 'G'}
	};

	std::map<char, int> values_map = {
		{'C', 0},
		{'A', 1},
		{'T', 2},
		{'G', 3}
	};

	std::string generateComplementSequence(const char *seq, unsigned int seq_length) {	
		std::string result;

		for(int i = 0; i < seq_length; i++) {
			result += complements_map[seq[i]];
		}


		return result;
	}

	unsigned int convertKmerToInteger(const char* temp_seg_ptr, int j, unsigned int k) {
		unsigned int result = 0;

		for(int i = 0; i < k; i++) {
			unsigned int temp = (values_map[(temp_seg_ptr + j)[i]]) << (2 * (k - i - 1));
			result |= temp;
		}

		return result;
	}

	void searchForMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>, MerTupleHasher> &minimizers_set,
				unsigned int k, unsigned int window_length, unsigned int i_start, std::function<bool(unsigned int)> i_predicate,
				std::function<const char*(const char*, unsigned int)> segment_ptr_extractor, 
				std::function<int(unsigned int, bool)> j_start_calculator, std::function<bool(unsigned int, int)> j_predicate, std::function<int(int)> j_mover,
				std::function<unsigned int(unsigned int, int)> min_pos_calculator) {
		
		long int position_of_last_found_min_mer = -1;
		unsigned int last_min_mer;

		for(unsigned int i = i_start; i_predicate(i); i++) {
			const char* temp_seg_ptr = segment_ptr_extractor(sequence, i);
			const char* temp_seg_complement_ptr = segment_ptr_extractor(sequence_complement, i);

			unsigned int min_mer = last_min_mer;
			bool is_min_complement;
			unsigned int min_pos = position_of_last_found_min_mer;
			bool first_found = position_of_last_found_min_mer >= i;

			for(int j = j_start_calculator(i, first_found); j_predicate(i, j); j = j_mover(j)) {

				unsigned int temp_mer_not_complement = convertKmerToInteger(temp_seg_ptr, j, k);
				unsigned int temp_mer_complement = convertKmerToInteger(temp_seg_complement_ptr, j, k);

				bool is_complement = temp_mer_complement < temp_mer_not_complement;
				unsigned int temp_mer = is_complement ? temp_mer_complement : temp_mer_not_complement;

				if(!first_found || min_mer > temp_mer) {
					if(!first_found) first_found = true;

					min_mer = temp_mer;
					is_min_complement = is_complement;
					min_pos = min_pos_calculator(i, j);
				}
			}

			minimizers_set.emplace(min_mer, min_pos, is_min_complement);
			position_of_last_found_min_mer = min_pos;
			last_min_mer = min_mer;
		}

	}

	void endMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>, MerTupleHasher> &minimizers_set, unsigned int k, unsigned int window_length, bool isStart) {
		if(isStart) {
			searchForMinimizers(sequence, sequence_complement, sequence_length, minimizers_set, k, window_length,
					k, [&](unsigned int i){ return i < k + window_length - 1;}, [](const char* seq, unsigned int i){ return seq;},
					[&](unsigned int i, bool condition) {return condition ? i - k : 0;}, [&](unsigned int i, int j){ return j + k <= i;},
					[&](int j){ return j + 1;}, [](unsigned int i, int j){ return j;});
		} else {
			searchForMinimizers(sequence, sequence_complement, sequence_length, minimizers_set, k, window_length,
					k, [&](unsigned int i){ return i < k + window_length - 1;}, [&](const char* seq, unsigned int i){ return seq + sequence_length - i;},
					[&](unsigned int i, bool condition) {return condition ? 0 : i - k;}, [](unsigned int i, int j){ return j >= 0;},
					[&](int j){ return j - 1;}, [&](unsigned int i, int j){ return sequence_length - i + j;});
		}		
	}

	void middleMinimizers(const char* sequence, const char* sequence_complement, unsigned int sequence_length, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>, MerTupleHasher> &minimizers_set, unsigned int k, unsigned int window_length) { 

		searchForMinimizers(sequence, sequence_complement, sequence_length, minimizers_set, k, window_length,
					0, [&](unsigned int i){ return window_length + k + i - 1 <= sequence_length;}, [](const char* seq, unsigned int start){ return seq + start;},
					[&](unsigned int i, bool condition) {return condition ? window_length - 1 : 0;}, [&](unsigned int i, int j){ return j + k <= window_length + k - 1;},
					[](int j){ return j + 1;}, [](unsigned int i, int j){ return i + j;});
	}

	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length) {
		
		std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vec;
		std::unordered_set<std::tuple<unsigned int, unsigned int, bool>, MerTupleHasher> minimizers_set;

		const char* seq_complement = generateComplementSequence(sequence, sequence_length).c_str();

		middleMinimizers(sequence, seq_complement, sequence_length, minimizers_set, k, window_length);

		endMinimizers(sequence, seq_complement, sequence_length, minimizers_set, k, window_length, true);

		endMinimizers(sequence, seq_complement, sequence_length, minimizers_set, k, window_length, false);

		minimizers_vec.assign(minimizers_set.begin(), minimizers_set.end());
		return minimizers_vec;
	}
}
