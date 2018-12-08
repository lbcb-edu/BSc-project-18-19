#include <map>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>

namespace orange {

	std::map<char, char> complements_map = {
		{A, T},
		{T, A},
		{G, C},
		{C, G}
	}

	std::string constructComplementMer(std::string const mer) {
		std:string complement;

		for(int i = 0; i < mer.length(); i++) {
			complement += complements_map(mer.at(i));
		}

		return complement;
	}

	std::map<unsigned int, std::pair<std::string, bool>> minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length) {
		std::string seq;
		seq.assign(sequence, sequence_length);
		
		std::map<unsigned int, std::pair<std::string, bool>> minimizers_map;
		unsigned int position_of_last_found_min_mer = 0;

		for(unsigned int i = 0; window_length + k + i - 1 <= seq.length(); i++) {
			std::string temp_seg = seq.substr(i, window_length + k + i - 1);

			std::string min_mer;
			bool is_min_complement;
			unsigned int min_pos = position_of_last_found_min_mer;
			bool first_found = false;

			for(int j = position_of_last_found_min_mer - i; j + k < temp_seg.length(); j++) {
				if(j < 0) j = 0;

				std::string temp_mer_not_complement = seq.substr(j, k);
				std::string temp_mer_complement = constructComplementMer(temp_mer);

				std::string temp_mer;
				bool is_complement;

				if(temp_mer_complement.compare(temp_mer_not_complement) < 0) {
					temp_mer = temp_mer_complement;
					is_complement = true;
				} else {
					temp_mer = temp_mer_not_complement;
					is_complement = false;
				}

				if(!first_found || min_mer.compare(temp_mer) > 0) {
					if(!first_found) first_found = true;

					min_mer = temp_mer;
					is_min_complement = is_complement;
					min_pos = i + j;
				}
			}

			if(minimizers_map.count(min_pos) == 0) {
				minimizers_map.insert(std::pair<unsigned int, std::pair<std::string, bool>>(min_pos, std::pair<std::string, bool>(min_mer, is_min_complement)));
				position_of_last_found_min_mer = min_pos;
			}
		}

		return minimizers_map;
	}
}
