#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <algorithm>

#include "white_minimizers.h"

namespace white {

char* reverseComplement(const char* sequence, unsigned int sequence_length) {
	char* newSequence = new char[sequence_length + 1];

	for (unsigned int i = 0; i < sequence_length; i++) {
		switch (sequence[i]) {
		case 'C':
			newSequence[sequence_length - 1 - i] = 'G';
			break;
		case 'A':
			newSequence[sequence_length - 1 - i] = 'T';
			break;
		case 'T':
			newSequence[sequence_length - 1 - i] = 'A';
			break;
		case 'G':
			newSequence[sequence_length - 1 - i] = 'C';
			break;
		default:
			std::cout << "Found a letter that doesn't belong to a DNA sequence.\n";
			exit(1);
		}
	}

	newSequence[sequence_length] = '\0';
	return newSequence;
}

std::vector <char*> create_kmers(const char* sequence, unsigned int sequence_length, unsigned int k) {
	std::vector <char*> k_mers;

	for (unsigned int i = 0; i <= sequence_length - k; i++) {
		char* newSequence = new char[k+1];
		for (unsigned int j = 0; j < k; j++) {
			newSequence[j] = sequence[i + j];
		}
		newSequence[k] = '\0';
		k_mers.push_back(newSequence);
	}

	return k_mers;
}

/*char* changeBasesToNumbers(const char* sequence, unsigned int sequence_length) {
	char* newSequence = new char[sequence_length + 1];

	for (unsigned int i = 0; i < sequence_length; i++) {
		if ((i+1) % 2 == 1)
			switch (sequence[i]) {
			case 'C':
				newSequence[i] = '1';
				break;
			case 'A':
				newSequence[i] = '2';
				break;
			case 'T':
				newSequence[i] = '3';
				break;
			case 'G':
				newSequence[i] = '4';
				break;
			default:
				std::cout << "Found a letter that doesn't belong to a DNA sequence.\n";
				exit(1);
			}
		else
			switch (sequence[i]) {
			case 'G':
				newSequence[i] = '1';
				break;
			case 'T':
				newSequence[i] = '2';
				break;
			case 'A':
				newSequence[i] = '3';
				break;
			case 'C':
				newSequence[i] = '4';
				break;
			default:
				std::cout << "Found a letter that doesn't belong to a DNA sequence.\n";
				exit(1);
			}
	}
	newSequence[sequence_length] = '\0';
	return newSequence;
}*/

char* changeBasesToNumbers(const char* sequence, unsigned int sequence_length) {
	char* newSequence = new char[sequence_length + 1];

	for (unsigned int i = 0; i < sequence_length; i++) {
			switch (sequence[i]) {
			case 'A':
				newSequence[i] = '1';
				break;
			case 'C':
				newSequence[i] = '2';
				break;
			case 'G':
				newSequence[i] = '3';
				break;
			case 'T':
				newSequence[i] = '4';
				break;
			default:
				std::cout << "Found a letter that doesn't belong to a DNA sequence.\n";
				exit(1);
			}

	}

	newSequence[sequence_length] = '\0';
	return newSequence;
}

int compare_kmers(char* kmer1, char* kmer2, unsigned int k) {
	int result = 0;

	for (unsigned int i = 0; i < k; i++) {
		if (kmer1[i] > kmer2[i]) {
			result = 1;
			break;
		}
		else if (kmer1[i] < kmer2[i]) {
			result = -1;
			break;
		}
	}
	return result;
}

std::tuple<unsigned int, unsigned int, bool> findMinimizerInWindow(
	std::vector <char*> window_from_original_sequence,
	std::vector <char*> window_from_reverse_complement,
	unsigned int starting_position,
	unsigned int k) {

	char* min = window_from_original_sequence.at(0);
	bool isOriginal = true;
	int counter_original = 0;
	int counter_reverse_complement = 0;

	for (unsigned int i = 1; i < window_from_original_sequence.size(); i++) {
		if (compare_kmers(min, window_from_original_sequence.at(i), k) == 1) {
			min = window_from_original_sequence.at(i);
			counter_original = i;
		}
	}

	for (unsigned int i = 0; i < window_from_reverse_complement.size(); i++) {
		if (compare_kmers(min, window_from_reverse_complement.at(i), k) == 1) {
			min = window_from_reverse_complement.at(i);
			counter_reverse_complement = i;
			isOriginal = false;
		}
	}

	return isOriginal ? std::make_tuple(atoi(min), starting_position + counter_original, true) :
		std::make_tuple(atoi(min), starting_position + counter_reverse_complement, false);
}

std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(
	const char* sequence,
	unsigned int sequence_length,
	unsigned int k,
	unsigned int window_length) {

	if (k > sequence_length) {
		std::cout << "K is bigger than the actual sequence." << std::endl;
		exit(1);
	}
	else if (window_length > sequence_length - k + 1) {
		std::cout << "Window size is bigger than the number of k-mers." << std::endl;
		exit(1);
	}

	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;
	std::tuple<unsigned int, unsigned int, bool> current_minimizer;

	std::vector <char*> current_window_original;
	std::vector <char*> current_window_reverse_complement;

	char* original_sequence = changeBasesToNumbers(sequence, sequence_length);
	char* temp = reverseComplement(sequence, sequence_length);
	char* reverse_complement = changeBasesToNumbers(temp, sequence_length);
	delete[] temp;

	std::vector <char*> k_mers_original = create_kmers(original_sequence, sequence_length, k);
	std::vector <char*> k_mers_reverse_complement = create_kmers(reverse_complement, sequence_length, k);

	for (unsigned int i = 0; i < window_length; i++) {
		current_window_original.push_back(k_mers_original.at(i));
		current_window_reverse_complement.push_back(k_mers_reverse_complement.at(i));
	}

	for (unsigned int i = 0; i < sequence_length - (window_length - 1 + k) + 1; i++) {
		current_minimizer = findMinimizerInWindow(current_window_original, current_window_reverse_complement, i, k);

		if (std::find(minimizers.begin(), minimizers.end(), current_minimizer) == minimizers.end()) {
			minimizers.push_back(current_minimizer);
		}

		if (window_length + i < k_mers_original.size()) {
			current_window_original.erase(current_window_original.begin());
			current_window_original.push_back(k_mers_original.at(window_length + i));

			current_window_reverse_complement.erase(current_window_reverse_complement.begin());
			current_window_reverse_complement.push_back(k_mers_reverse_complement.at(window_length + i));
		}

	}

	delete[] original_sequence;
	delete[] reverse_complement;

	for (auto k_mer : k_mers_original) {
		delete[] k_mer;
	}
	for (auto k_mer : k_mers_reverse_complement) {
		delete[] k_mer;
	}

	return minimizers;
}

/*int main()
{
	std::string seq = std::string("GTCATGCACGTTCAC");
	const char* sequence = seq.c_str();

	std::vector<std::tuple<unsigned int, unsigned int, bool>> mins = minimizers(sequence, 15, 3, 4);

	return 0;
}*/

}
