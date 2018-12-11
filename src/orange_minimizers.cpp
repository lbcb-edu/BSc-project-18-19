#include <map>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <tuple>
#include <algorithm>
#include <functional>

namespace orange {

	std::map<char, char> complements_map = {
		{'A', 'T'},
		{'T', 'A'},
		{'G', 'C'},
		{'C', 'G'}
	};

	std::map<char, char> values_map = {
		{'C', '0'},
		{'A', '1'},
		{'T', '2'},
		{'G', '3'}
	};

	
	std::string constructStringWithMap(std::string const base_string, std::map<char, char> const map) {
		std::string result;

		for(int i = 0; i < base_string.length(); i++) {
			result += map.at(base_string.at(i));
		}

		return result;
	}

	std::string constructComplementMer(std::string const mer) {
		return constructStringWithMap(mer, complements_map);
	}

	std::string constructValueString(std::string const mer) {
		return constructStringWithMap(mer, values_map);
	}

	void searchForMinimizers(const std::string seq, std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length,
				unsigned int i_start, std::function<bool(unsigned int)> i_predicate, std::function<std::string(unsigned int)> segment_extractor,
				std::function<int(unsigned int, bool)> j_start_calculator, std::function<bool(unsigned int, int)> j_predicate, std::function<int(int)> j_mover,
				std::function<unsigned int(unsigned int, int)> min_pos_calculator) {
		
		unsigned int position_of_last_found_min_mer = -1;
		unsigned int last_min_mer;

		for(unsigned int i = i_start; i_predicate(i); i++) {
			std::string temp_seg = segment_extractor(i);

			unsigned int min_mer = last_min_mer;
			bool is_min_complement;
			unsigned int min_pos = position_of_last_found_min_mer;
			bool first_found = (position_of_last_found_min_mer >= i && position_of_last_found_min_mer <= window_length + k + i - 1);

			for(int j = j_start_calculator(i, first_found); j_predicate(i, j); j = j_mover(j)) {

				unsigned int temp_mer_not_complement = strtol(constructValueString(temp_seg.substr(j, k)).c_str(), NULL, 4);
				unsigned int temp_mer_complement = strtol(constructValueString(constructComplementMer(temp_seg.substr(j, k))).c_str(), NULL, 4);

				bool is_complement = temp_mer_complement < temp_mer_not_complement;
				unsigned int temp_mer = is_complement ? temp_mer_complement : temp_mer_not_complement;

				if(!first_found || min_mer > temp_mer) {
					if(!first_found) first_found = true;

					min_mer = temp_mer;
					is_min_complement = is_complement;
					min_pos = min_pos_calculator(i, j);
				}
			}

			std::tuple<unsigned int, unsigned int, bool> temp_tuple = std::make_tuple (min_mer, min_pos, is_min_complement);

			if(std::find(minimizers_vec.begin(), minimizers_vec.end(), temp_tuple) == minimizers_vec.end()) {
				minimizers_vec.push_back(temp_tuple);
				position_of_last_found_min_mer = min_pos;
				last_min_mer = min_mer;
			}
		}
	}

	void endMinimizers(const std::string seq, std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length, bool isStart) {
		if(isStart) {
			searchForMinimizers(seq, minimizers_vec, k, window_length,
					k, [&](unsigned int i){ return i < k + window_length - 1;}, [&](unsigned int i){ return seq.substr(0, i);},
					[&](unsigned int i, bool condition) {return condition ? i - k : 0;}, [&](unsigned int i, int j){ return j + k <= i;},
					[&](int j){ return j + 1;}, [](unsigned int i, int j){ return j;});
		} else {
			searchForMinimizers(seq, minimizers_vec, k, window_length,
					k, [&](unsigned int i){ return i < k + window_length - 1;}, [&](unsigned int i){ return seq.substr(seq.length() - i, i);},
					[&](unsigned int i, bool condition) {return condition ? 0 : i - k;}, [](unsigned int i, int j){ return j >= 0;},
					[&](int j){ return j - 1;}, [&](unsigned int i, int j){ return seq.length() - i + j;});
		}		
	}

	void middleMinimizers(const std::string seq, std::vector<std::tuple<unsigned int, unsigned int, bool>> &minimizers_vec, unsigned int k, unsigned int window_length) { 

		searchForMinimizers(seq, minimizers_vec, k, window_length,
					0, [&](unsigned int i){ return window_length + k + i - 1 <= seq.length();}, [&](unsigned int start){ return seq.substr(start, window_length + k - 1);},
					[&](unsigned int i, bool condition) {return condition ? window_length - 1 : 0;}, [&](unsigned int i, int j){ return j + k <= window_length + k - 1;},
					[](int j){ return j + 1;}, [](unsigned int i, int j){ return i + j;});
	}

	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length) {
		std::string seq;
		seq.assign(sequence, sequence_length);
		
		std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vec;

		middleMinimizers(seq, minimizers_vec, k, window_length);

		endMinimizers(seq, minimizers_vec, k, window_length, true);

		endMinimizers(seq, minimizers_vec, k, window_length, false);

		return minimizers_vec;
	}
}
