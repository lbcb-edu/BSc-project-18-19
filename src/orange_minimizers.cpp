#include <map>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>

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

	std::map<unsigned int, std::pair<std::string, bool>> minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length) {
		std::string seq;
		seq.assign(sequence, sequence_length);
		
		std::map<unsigned int, std::pair<std::string, bool>> minimizers_map;
		unsigned int position_of_last_found_min_mer = -1;
		std::string last_min_mer;

		for(unsigned int i = 0; window_length + k + i - 1 <= seq.length(); i++) {
			std::string temp_seg = seq.substr(i, window_length + k - 1);

			std::string min_mer = last_min_mer;
			bool is_min_complement;
			unsigned int min_pos = position_of_last_found_min_mer;
			bool first_found = (position_of_last_found_min_mer >= i && position_of_last_found_min_mer <= window_length + k + i - 1);

			for(int j = first_found ? temp_seg.length() - k : 0; j + k <= temp_seg.length(); j++) {

				std::string temp_mer_not_complement = temp_seg.substr(j, k);
				std::string temp_mer_complement = constructComplementMer(temp_mer_not_complement);

				bool is_complement = constructValueString(temp_mer_complement).compare(constructValueString(temp_mer_not_complement)) < 0;
				std::string temp_mer = is_complement ? temp_mer_complement : temp_mer_not_complement;

				if(!first_found || constructValueString(min_mer).compare(constructValueString(temp_mer)) > 0) {
					if(!first_found) first_found = true;

					min_mer = temp_mer;
					is_min_complement = is_complement;
					min_pos = i + j;
				}
			}

			if(minimizers_map.count(min_pos) == 0) {
				minimizers_map.insert(std::pair<unsigned int, std::pair<std::string, bool>>(min_pos, std::pair<std::string, bool>(min_mer, is_min_complement)));
				position_of_last_found_min_mer = min_pos;
				last_min_mer = min_mer;
			}
		}

		return minimizers_map;
	}
}
