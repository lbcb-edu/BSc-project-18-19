#include <map>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <tuple>
#include <algorithm>

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

	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length) {
		std::string seq;
		seq.assign(sequence, sequence_length);
		
		std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vec;
		unsigned int position_of_last_found_min_mer = -1;
		unsigned int last_min_mer;

		for(unsigned int i = 0; window_length + k + i - 1 <= seq.length(); i++) {
			std::string temp_seg = seq.substr(i, window_length + k - 1);

			unsigned int min_mer = last_min_mer;
			bool is_min_complement;
			unsigned int min_pos = position_of_last_found_min_mer;
			bool first_found = (position_of_last_found_min_mer >= i && position_of_last_found_min_mer <= window_length + k + i - 1);

			for(int j = first_found ? temp_seg.length() - k : 0; j + k <= temp_seg.length(); j++) {

				unsigned int temp_mer_not_complement = strtol(constructValueString(temp_seg.substr(j, k)).c_str(), NULL, 4);
				unsigned int temp_mer_complement = strtol(constructValueString(constructComplementMer(temp_seg.substr(j, k))).c_str(), NULL, 4);

				bool is_complement = temp_mer_complement < temp_mer_not_complement;
				unsigned int temp_mer = is_complement ? temp_mer_complement : temp_mer_not_complement;

				if(!first_found || min_mer > temp_mer) {
					if(!first_found) first_found = true;

					min_mer = temp_mer;
					is_min_complement = is_complement;
					min_pos = i + j;
				}
			}

			std::tuple<unsigned int, unsigned int, bool> temp_tuple = std::make_tuple (min_mer, min_pos, is_min_complement);

			if(std::find(minimizers_vec.begin(), minimizers_vec.end(), temp_tuple) == minimizers_vec.end()) {
				minimizers_vec.push_back(temp_tuple);
				position_of_last_found_min_mer = min_pos;
				last_min_mer = min_mer;
			}
		}

		return minimizers_vec;
	}
}
