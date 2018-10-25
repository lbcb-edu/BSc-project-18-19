#include<iostream>
#include<memory>
#include<vector>
#include<getopt.h>
#include<stdlib.h>
#include<set>
#include<string>

#include "white_mapper.h"
#include "bioparser/bioparser.hpp"

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
			const char* quality, uint32_t quality_length)
                        {
                                this -> name = std::string (name, name_length);
                                this -> sequence = std::string (sequence, sequence_length);
				this -> quality = std::string (quality, quality_length);
			}

		SequenceFormat (
                        const char* name, uint32_t name_length,
                        const char* sequence, uint32_t sequence_length
                        )
                        {
                                this -> name = std::string (name, name_length);
                                this -> sequence = std::string (sequence, sequence_length);
                                this -> quality = std::string ("");
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
				"\n	This program is a mapper that requires two files. First one containing a set		\n"
                     		"	of fragments in FASTA or FASTQ format and the second one containing corresponding	\n"
                                "	reference genome in FASTA format. The mapper parses the files and stores them in	\n"
                                "	memory while also providing statistics of given files:					\n"
                                "	-> number of sequences									\n"
                                "	-> average length									\n"
                                "	-> minimal and maximal length.								\n\n"

                                "	File extensions accepted:								\n"
                                "	-> .fasta             -> .fastq							\n"
                                "	-> .fa                -> .fq								\n"
                                "	-> .fasta.gz          -> .fastq.gz							\n"
                                "	-> .fa.gz             -> .fq.gz							\n\n"

				"	Usage: white_mapper [OPTIONS] sequences_file reference_genome_file			\n\n"

                                "	There are two options:									\n"
                                "	-> \"-h\" or \"--help\" for help							\n"
                                "	-> \"-v\" or \"--version\" for displaying the current version.				\n\n";
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
void parse_file (std::string file_path)
{
	std::vector<std::unique_ptr<SequenceFormat>> objects;
	auto parser = bioparser::createParser<T, SequenceFormat> (file_path);
	parser -> parse_objects (objects, -1);

	statistics (objects, file_path);
}

int main (int argc, char* argv[])
{
	int option = 0;
	int option_done[] = {0, 0};

        while ((option = getopt_long(argc, argv, "hv", long_options, NULL)) != -1) {
                switch (option) {
                        case 'h':

				if(option_done[0] == 0) {
					option_done[0] = 1;
                                	help();
				}
                                break;

                        case 'v':

				if(option_done[1] == 0) {
                                        option_done[1] = 1;
                                        version();
                                }
                                break;

                        default:
                                fprintf (stderr, "\n--Error:\n\tUnsupported option. Usage: %s [-h] [--help] [-v] [--version] sequence_file reference_genome_file\n\n", argv[0]);
                                break;


                }
	}

	if (optind == argc && option != 0)
		exit (0);


	if (argc - optind == 2)
	{
		if (check_format (argv[optind], fasta_formats))
			parse_file<bioparser::FastaParser> (std::string (argv[optind]));

		else if (check_format (argv[optind], fastq_formats))
                        parse_file<bioparser::FastqParser> (std::string (argv[optind]));

		else
			fprintf (stderr, "\nError:													\n"
					 "	Unsuported format. File containing sequences (1st file) needs to have one of the following extensions:	\n"
                                     	 "	-> .fasta             -> .fastq										\n"
                                     	 "	-> .fa                -> .fq										\n"
                                     	 "	-> .fasta.gz          -> .fastq.gz									\n"
                                     	 "	-> .fa.gz             -> .fq.gz										\n\n");

		if (check_format (argv[++optind], fasta_formats))
                        parse_file<bioparser::FastaParser> (std::string (argv[optind]));

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
					 "	Usage: white_mapper [OPTIONS] sequence_file reference_genome_file	\n\n"

					 "	Number of nonoption arguments must be 2:				\n"
                                         "		1st file contains sequences in FASTA or FASTQ format		\n"
                                         "		2nd file contains a reference genome in FASTA fromat		\n\n");

	return 0;
}
