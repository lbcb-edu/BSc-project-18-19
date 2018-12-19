#include <iostream>
#include <memory>
#include <vector>
#include <getopt.h>
#include <stdlib.h>
#include <set>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <cmath>

#include "white_minimizers.h"
#include "white_alignment.h"
#include "white_mapper.h"
#include "bioparser/bioparser.hpp"

#define MATCH 2;
#define MISMATCH -1;
#define GAP -2;
#define K_MER_LENGTH 15;
#define WINDOW_SIZE 5;
#define PERCENTAGE 0.001;

class SequenceFormat
{
	 public:
                std::string name;
                std::string sequence;
		std::string quality;

		friend bioparser::FastaParser<SequenceFormat>;
		friend bioparser::FastqParser<SequenceFormat>;
		friend void statistics (std::vector<std::unique_ptr<SequenceFormat>>&);

	private:
		SequenceFormat (
                        const char* name, uint32_t name_length,
                        const char* sequence, uint32_t sequence_length,
			const char* quality, uint32_t quality_length
			) : name (name, name_length), sequence (sequence, sequence_length), quality (quality, quality_length)
                        {
			}

		SequenceFormat (
                        const char* name, uint32_t name_length,
                        const char* sequence, uint32_t sequence_length
                        ) : name (name, name_length), sequence (sequence, sequence_length), quality ("")
                        {
                        }

};

const std::set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};


const struct option long_options[] = {
{"help",	no_argument,	0,	'h' },
{"version",	no_argument,	0,	'v' },
{0,		0,		0,	 0  },
};

void statistics (std::vector<std::unique_ptr<SequenceFormat>> &v, std::string file)
{
        uint32_t min_length = v.at(0) -> sequence.size();
        uint32_t max_length = min_length;
        uint32_t total_length = 0;

        for (auto &ptr : v)
        {
                if (ptr -> sequence.size() < min_length)
                        min_length = ptr -> sequence.size();

                if (ptr -> sequence.size() > max_length)
                        max_length = ptr -> sequence.size();

                total_length += ptr -> sequence.size();
        }

	printf( "\nStatistics for file: %s\n"
		"	Number of sequences: %ld\n"
		"	Minimum length: %d\n"
		"	Maximum length: %d\n"
		"	Average length: %.2f\n\n" ,
		file.c_str(), v.size(), min_length, max_length, (float)total_length/v.size() );

}

void help() {
	std::cout << 		"--Help:"
				"\n	This program is a mapper that requires two files. First file should contain a set of	\n"
				"	fragments in FASTA or FASTQ format and the second one a corresponding reference genome	\n"
                                "	in FASTA format. The mapper parses the files and stores them in	memory while also	\n"
				"	providing statistics of given files:							\n"
                                "	-> number of sequences									\n"
                                "	-> average length									\n"
                                "	-> minimal and maximal length.								\n\n"

				"	Afterwards, two random sequences from the first file are selected and aligned using	\n"
				"	global, semi-global and local alignment. The resulting cigar string is printed.		\n"
				"	Default for match, mismatch ang gap values is 2, -1, -2 respectively. The default	\n"
				"	values can be changed using appropriate options.					\n\n"

				"	Seed and extend:									\n"
				"	In order to alleviate the alignment proces we use seed and extend approach. Among all	\n"
				"	k-mers we choose a set of minimizers which will be used for fast detection of similar	\n"
				"	regions between two sequences prior the exact alignment. The mapper will find minimizers\n"
				"	for every sequence in the first file and print the number of distinct minimizers and the\n"
				"	number of occurences of the most frequent minimizer when the top f frequent minimizers 	\n"
				"	are not taken in account.								\n\n"

				"	The proces of finding minimizers uses three variables:					\n"
				"	-> k-mer length k									\n"
				"	-> window size w									\n"
				"	-> percentage of top minimizers to disregard f						\n\n"

				"	Their default values are 15, 5 and 0.001 respecively.					\n\n"

                                "	File extensions accepted:								\n"
                                "	-> .fasta             -> .fastq								\n"
                                "	-> .fa                -> .fq								\n"
                                "	-> .fasta.gz          -> .fastq.gz							\n"
                                "	-> .fa.gz             -> .fq.gz								\n\n"

				"	Usage: white_mapper [OPTIONS] sequences_file reference_genome_file match mismatch gap	\n\n"

                                "	There are five options:									\n"
                                "	-> \"-h\" or \"--help\" for help							\n"
                                "	-> \"-v\" or \"--version\" for displaying the current version.				\n"
				"	-> \"-m ARG\" sets match value to ARG							\n"
				"       -> \"-s ARG\" sets mismatch value to ARG                                                \n"
				"       -> \"-g ARG\" sets gap value to ARG                                                     \n"
				"       -> \"-k ARG\" sets length k of k-mers to ARG						\n"
				"       -> \"-w ARG\" sets window size to ARG							\n"
				"	-> \"-f ARG\" sets f to ARG								\n\n";
}

void version() {
	std::cout << "\n--Version:";
	std::cout << " " <<white_mapper_VERSION_MAJOR <<"." <<white_mapper_VERSION_MINOR <<"\n\n";
}

bool check_format (std::string file, std::set<std::string> list_of_formats) {
	for (auto &format : list_of_formats)
		if(file.size() > format.size())
			if (file.compare (file.size() - format.size(), format.size(), format) == 0)
				return true;
	return false;
}

template < template<class> class T>
std::vector<std::unique_ptr<SequenceFormat>> parse_file (std::string file_path)
{
	std::vector<std::unique_ptr<SequenceFormat>> objects;
	auto parser = bioparser::createParser<T, SequenceFormat> (file_path);
	parser -> parse_objects (objects, -1);

	statistics (objects, file_path);

	return objects;
}

void minimizer_occurrences (std::vector<std::unique_ptr<SequenceFormat>> &sequences, unsigned int k, unsigned int window_size, float f) {
	std::map <unsigned int, unsigned int> minimizer_occurrences;
	std::vector <std::tuple<unsigned int, unsigned int, bool>> current_minimizers;

	for (auto &ptr : sequences) {
		current_minimizers = white::minimizers(ptr -> sequence.c_str(), ptr -> sequence.size(), k, window_size);

		for (auto minimizer_tuple : current_minimizers) {
			minimizer_occurrences[std::get<0>(minimizer_tuple)]++;
		}
	}

	/*std::ofstream fout;
  	fout.open ("minimizer_occurrences.csv", std::ios::out);

	if (!fout.is_open()){
		std::cout << "Unable to open file.\n";
		exit(1);
	}

  	fout << "Minimizer,Number of occurrences\r\n"; //posto mi je linux na windowsima da mogu citat i u notepadu

	for (std::map<unsigned int, unsigned int>::iterator it = minimizer_occurrences.begin(); it != minimizer_occurrences.end(); it++)
	{
		if(it->second != 1) {
  			fout << it->first << ",";
			fout << it->second << "\r\n";
		}
	}

  	fout.close();*/

	unsigned int number_of_minimizers_to_disregard = (unsigned int) std::round(minimizer_occurrences.size() * f);

	std::map<unsigned int, unsigned int>::iterator it = std::prev(minimizer_occurrences.end(), number_of_minimizers_to_disregard + 1);

	printf("\n--- Minimizer statistics ---\n\n"

		"-> Minimizers found: %lu								\n"
		"-> After disregarding top %f minimizers the most frequent minimizer is:		\n"
		"	Minimizer: %u									\n"
		"	Frequency: %u									\n\n",

		minimizer_occurrences.size(), f, it->first, it->second);
}

int main (int argc, char* argv[])
{
	int option = 0;
	std::vector<std::unique_ptr<SequenceFormat>> sequences;
	std::vector<std::unique_ptr<SequenceFormat>> reference_genome;

	int match = MATCH;
        int mismatch = MISMATCH;
        int gap = GAP;
	int k = K_MER_LENGTH;
	int window_size = WINDOW_SIZE;
	float f = PERCENTAGE;

        while ((option = getopt_long(argc, argv, "hvm:s:g:k:w:f:", long_options, NULL)) != -1) {
                switch (option) {
                        case 'h':
				help();
                                exit(0);

                        case 'v':
				version();
                                exit(0);

			case 'm':
				match = atoi (optarg);
				break;

			case 's':
                                mismatch = atoi (optarg);
				break;

			case 'g':
                                gap = atoi (optarg);
				break;

			case 'k':
				k = atoi (optarg);
				break;

			case 'w':
				window_size = atoi (optarg);
				break;

			case 'f':
				f = std::stof (optarg, NULL);
				break;

                        default:
                                fprintf (stderr, "\n--Error:\n\tUnsupported option. Usage: %s [-h] [--help] [-v] [--version] [-m ARG] [-s ARG] [-g ARG] [-k ARG] [-w ARG] [-f ARG] "
						"sequence_file reference_genome_file match mismatch gap\n\n", argv[0]);
                                exit(0);

                }
	}

	if (argc - optind == 2)
	{
		if (check_format (argv[optind], fasta_formats))
			sequences = parse_file<bioparser::FastaParser> (std::string (argv[optind]));

		else if (check_format (argv[optind], fastq_formats))
                        sequences = parse_file<bioparser::FastqParser> (std::string (argv[optind]));

		else
			fprintf (stderr, "\nError:													\n"
					 "	Unsuported format. File containing sequences (1st file) needs to have one of the following extensions:	\n"
                                     	 "	-> .fasta             -> .fastq										\n"
                                     	 "	-> .fa                -> .fq										\n"
                                     	 "	-> .fasta.gz          -> .fastq.gz									\n"
                                     	 "	-> .fa.gz             -> .fq.gz										\n\n");

		if (check_format (argv[++optind], fasta_formats))
                        reference_genome = parse_file<bioparser::FastaParser> (std::string (argv[optind]));
		else

			fprintf (stderr, "\n --Error:														\n"
					 "	Unsuported format. File containing reference genome (2nd file) needs to have one of the following extensions:	\n"
                                     	 "	-> .fasta													\n"
                                     	 "	-> .fa														\n"
                                     	 "	-> .fasta.gz													\n"
                                     	 "	-> .fa.gz													\n\n");
	}

	 else
                        fprintf (stderr, "\n--Error:\n"
					 "	Usage: white_mapper [OPTIONS] sequence_file reference_genome_file match mismatch gap	\n\n"

					 "	Number of nonoption arguments must be 2:						\n"
                                         "		1st file contains sequences in FASTA or FASTQ format				\n"
                                         "		2nd file contains a reference genome in FASTA fromat				\n\n");

	//to use these next few lines change constructor of SequenceFormat to public and minimizer_occurrences (sekvenceTest, k, window_size)
	/*std::vector<std::unique_ptr<SequenceFormat>> sekvenceTest;
	std::unique_ptr<SequenceFormat> p1 (new SequenceFormat ("S1", 2, "TCAGGAAGAAGCAGA", 15));
	std::unique_ptr<SequenceFormat> p2 (new SequenceFormat ("S2", 2, "GTCATGCACGTTCAC", 15));
	std::unique_ptr<SequenceFormat> p3 (new SequenceFormat ("S3", 2, "TCAGGAAGAAGCAGA", 15));
	sekvenceTest.push_back(std::move(p1));
	sekvenceTest.push_back(std::move(p2));
	sekvenceTest.push_back(std::move(p3));*/

	//finding minimizers and making a csv file of their occurrences
	minimizer_occurrences (sequences, k, window_size, f);

	//alignment of 2 random selected sequences
	srand (time (NULL));

	int i1 = rand() % sequences.size();
	int i2 = rand() % sequences.size();

	//printf("Sequences used:\n%d_%d", i1, i2);

	int values[3];

	std::string cigar;
	unsigned int target_begin = 0;
//	printf("\nFirst sequence (Length: %lu):\n%s\n\nSecond sequence (Length: %lu):\n%s\n\n", sequences[i1] -> sequence.size(), sequences[i1] -> sequence.c_str(),
//						sequences[i2] -> sequence.size(), sequences[i2] -> sequence.c_str());

	int value = white::pairwise_alignment(	sequences[i1] -> sequence.c_str(), sequences[i1] -> sequence.size(),
  	                                   	sequences[i2] -> sequence.c_str(), sequences[i2] -> sequence.size(),
 	                                  	white::AlignmentType::global, match, mismatch, gap, cigar, target_begin);

	values[0] = value;

	std::cout << "\nGlobal alignment:\n";
	//std::cout << value << std::endl;
	//std::cout << target_begin << std::endl;
	printf ("%s\n", cigar.c_str());

	value = white::pairwise_alignment(  sequences[i1] -> sequence.c_str(), sequences[i1] -> sequence.size(),
                                            sequences[i2] -> sequence.c_str(), sequences[i2] -> sequence.size(),
                                            white::AlignmentType::semi_global, match, mismatch, gap, cigar, target_begin);

	values[1] = value;

        std::cout << "\nSemi-global alignment:\n";
        //std::cout << value << std::endl;
        //std::cout << target_begin << std::endl;
	printf ("%s\n", cigar.c_str());

	value = white::pairwise_alignment(  sequences[i1] -> sequence.c_str(), sequences[i1] -> sequence.size(),
                                            sequences[i2] -> sequence.c_str(), sequences[i2] -> sequence.size(),
                                            white::AlignmentType::local, match, mismatch, gap, cigar, target_begin);

	values[2] = value;

        std::cout << "\nLocal alignment:\n";
        //std::cout << value << std::endl;
        //std::cout << target_begin << std::endl;
	printf ("%s\n", cigar.c_str());


/*	std::string q = {"TCCG"};
	std::string t = {"ACTCCGAT"};
	std::string cigar;
	unsigned int target_begin = 0;
	int value = white::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), white::AlignmentType::semi_global, match, mismatch, gap, cigar, target_begin);
	std::cout << value << std::endl;
	std::cout << target_begin << std::endl;
	std::cout << cigar << std::endl;
	std::cout <<std::endl;
*/
/*	std::string q = {"TTCCGCCAA"};
	std::string t = {"AACCCCTT"};
	std::string cigar;
	unsigned int target_begin = 0;
	int value = white::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), white::AlignmentType::local, match, mismatch, gap, cigar, target_begin);
	std::cout << value << std::endl;
	std::cout << target_begin << std::endl;
	std::cout << cigar << std::endl;
	std::cout <<std::endl;
*/

/*	std::string q = {"TGCATAT"};
	std::string t = {"ATCCGAT"};
	std::string cigar;
	unsigned int target_begin = 0;
	int value = white::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), white::AlignmentType::global, match, mismatch, gap, cigar, target_begin);
	std::cout << value << std::endl;
	std::cout << target_begin << std::endl;
	std::cout << cigar << std::endl;
	std::cout <<std::endl;
*/
	return 0;
}
