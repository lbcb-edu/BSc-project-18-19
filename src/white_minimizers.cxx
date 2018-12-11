#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <algorithm>
#include <map>
#include <fstream>
#include <memory>

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

	unsigned int my_atoi(const char* string, unsigned int number_of_chars) {
		int result = 0;

		for (unsigned int i = 0; i < number_of_chars; i++)
			result = result * 10 + string[i] - '0';

		return result;
	}

	int compare_kmers(const char* kmer1, const char* kmer2, unsigned int k) {
		unsigned int k1 = my_atoi(kmer1, k);
		unsigned int k2 = my_atoi(kmer2, k);

		if (k1 > k2)
			return 1;
		else if (k1 < k2)
			return -1;
		else
			return 0;
	}

	std::tuple<unsigned int, unsigned int, bool> findMinimizerInWindow(
		const char* original_sequence,
		const char* reverse_complement,
		unsigned int starting_position,
		unsigned int k,
		unsigned int window_size) {

		unsigned int min_original = starting_position;
		unsigned int min_reverse_complement = starting_position;

		for (unsigned int i = 1; i < window_size; i++) {
			if (compare_kmers(&original_sequence[min_original], &original_sequence[starting_position + i], k) == 1)
				min_original = starting_position + i;
		}

		for (unsigned int i = 1; i < window_size; i++) {
			if (compare_kmers(&reverse_complement[min_reverse_complement], &reverse_complement[starting_position + i], k) == 1)
				min_reverse_complement = starting_position + i;
		}

		return compare_kmers(&original_sequence[min_original], &reverse_complement[min_reverse_complement], k) <= 0 ?
			std::make_tuple(my_atoi(&original_sequence[min_original], k), min_original, true) :
			std::make_tuple(my_atoi(&reverse_complement[min_reverse_complement], k), min_reverse_complement, false);
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

		char* original_sequence = changeBasesToNumbers(sequence, sequence_length);
		char* temp = reverseComplement(sequence, sequence_length);
		char* reverse_complement = changeBasesToNumbers(temp, sequence_length);
		delete[] temp;

		for (unsigned int i = 0; i < sequence_length - (window_length - 1 + k) + 1; i++) {
			current_minimizer = findMinimizerInWindow(original_sequence, reverse_complement, i, k, window_length);

			if (std::find(minimizers.begin(), minimizers.end(), current_minimizer) == minimizers.end()) {
				minimizers.push_back(current_minimizer);
			}
		}

		delete[] original_sequence;
		delete[] reverse_complement;

		return minimizers;
	}

	class SequenceFormat
	{
	public:
		std::string name;
		std::string sequence;
		std::string quality;

		SequenceFormat(
			const char* name, uint32_t name_length,
			const char* sequence, uint32_t sequence_length,
			const char* quality, uint32_t quality_length
		) : name(name, name_length), sequence(sequence, sequence_length), quality(quality, quality_length)
		{
		}

		SequenceFormat(
			const char* name, uint32_t name_length,
			const char* sequence, uint32_t sequence_length
		) : name(name, name_length), sequence(sequence, sequence_length), quality("")
		{
		}

	};

	void minimizer_occurrences(std::vector<std::unique_ptr<SequenceFormat>> &sequences, unsigned int k, unsigned int window_size) {
		std::map <unsigned int, unsigned int> minimizer_occurrences;
		std::vector <std::tuple<unsigned int, unsigned int, bool>> current_minimizers;

		for (auto &ptr : sequences) {
			current_minimizers = minimizers(ptr->sequence.c_str(), ptr->sequence.size(), k, window_size);

			for (auto minimizer_tuple : current_minimizers) {
				minimizer_occurrences[std::get<0>(minimizer_tuple)]++;
			}
		}

		std::ofstream fout;
		fout.open("minimizer_occurrences.csv");

		if (!fout.is_open()) {
			std::cout << "Unable to open file.\n";
			exit(1);
		}

		fout << "Minimizer, Number of occurrences\n";

		for (std::map<unsigned int, unsigned int>::iterator it = minimizer_occurrences.begin(); it != minimizer_occurrences.end(); it++)
		{
			fout << it->first << ", ";
			fout << it->second << "\n";
		}

		fout.close();
	}

	/*int main()
	{
		std::vector<std::unique_ptr<SequenceFormat>> sekvenceTest;
		std::unique_ptr<SequenceFormat> p1(new SequenceFormat("S1", 2, "TCAGGAAGAAGCAGA", 15));
		std::unique_ptr<SequenceFormat> p2(new SequenceFormat("S2", 2, "GTCATGCACGTTCAC", 15));
		std::unique_ptr<SequenceFormat> p3(new SequenceFormat("S3", 2, "TCAGGAAGAAGCAGA", 15));
		sekvenceTest.push_back(std::move(p1));
		sekvenceTest.push_back(std::move(p2));
		sekvenceTest.push_back(std::move(p3));
		//finding minimizers and making a csv file of their occurrences
		minimizer_occurrences(sekvenceTest, 3, 4);
		return 0;
	}*/

}
