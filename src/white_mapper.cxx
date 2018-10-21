#include<iostream>
#include<memory>
#include<vector>

#include "white_mapper.h"
#include "vendor/bioparser/include/bioparser/bioparser.hpp"

class SequenceFormat
{
	 public:
                const char* name;
                const char* sequence;
		const char* quality;
                uint32_t name_length;
                uint32_t sequence_length;
		uint32_t quality_length;
		friend bioparser::FastaParser<SequenceFormat>;
		friend bioparser::FastqParser<SequenceFormat>;
		friend void statistics (std::vector<std::unique_ptr<SequenceFormat>>&);

	private:
		SequenceFormat (
                        const char* name, uint32_t name_length,
                        const char* sequence, uint32_t sequence_length,
			const char* quality, uint32_t quality_length)
                        {
                                this -> name = name;
                                this -> name_length = name_length;
                                this -> sequence = sequence;
                                this -> sequence_length = sequence_length;
				this -> quality;
				this -> quality_length;
			}

		SequenceFormat (
                        const char* name, uint32_t name_length,
                        const char* sequence, uint32_t sequence_length
                        )
                        {
                                this -> name = name;
                                this -> name_length = name_length;
                                this -> sequence = sequence;
                                this -> sequence_length = sequence_length;
                                this -> quality = NULL;
                                this -> quality_length = 0;
                        }

};

void statistics (std::vector<std::unique_ptr<SequenceFormat>> &v)
{
        int min_length = v.at(0) -> sequence_length;
        int max_length = min_length;
        int total_length = 0;

        for (auto &ptr : v)
        {
                if (ptr -> sequence_length < min_length)
                        min_length = ptr -> sequence_length;

                if (ptr -> sequence_length > max_length)
                        max_length = ptr -> sequence_length;

                total_length += ptr -> sequence_length;
        }

	std::cout << "\nNumber of sequences: " << v.size() <<"\n";
        std::cout << "Minimum sequence length: " << min_length << "\n";
        std::cout << "Maximum sequence length: " << max_length << "\n";
        std::cout << "Average sequence length: " << (float) total_length / v.size() << "\n\n";

}


int main (int argc, char* argv[])
{
	if (argc == 2)
	{
		std::string option (argv[1]);

		if (option.compare("-v") == 0)
			std::cout << "\n\"" <<white_mapper_VERSION_MAJOR <<"." <<white_mapper_VERSION_MINOR <<"\"\n\n";

		else if (option.compare("-h") == 0)
			std::cout << "\nThis program is a mapper that requires two files. First one containing a set\n"
				"of fragments in FASTA or FASTQ format and the second one containing corresponding  \n"
				"reference genome in FASTA format. The mapper parses the files and stores them in   \n"
				"memory while also providing statistics of given files:                             \n"
				"-> number of sequences                                                             \n"
				"-> average length                                                                  \n"
				"-> minimal and maximal length.                                                   \n\n"

				"File extensions accepted:                                                          \n"
				"-> .fasta             -> .fastq                                                    \n"
				"-> .fa                -> .fq                                                       \n"
				"-> .fasta.gz          -> .fastq.gz                                                 \n"
				"-> .fa.gz             -> .fq.gz                                                  \n\n"

				"There are two options:                                                             \n"
				"-> \"-h\" for help                                                                 \n"
				"-> \"-v\" for displaying the current version.                                    \n\n";
		else
			std::cerr << "\nInvalid option \"" << argv[1] << "\".\n\n";
	}

	else if (argc == 3)
	{

		std::string sequences_file_path (argv[1]);
		std::string genome_file_path (argv[2]);

		if (
			sequences_file_path.rfind(".fasta") == sequences_file_path.size() - 6 ||
			sequences_file_path.rfind(".fa") == sequences_file_path.size() - 3 ||
			sequences_file_path.rfind(".fasta.gz") == sequences_file_path.size() - 9 ||
			sequences_file_path.rfind(".fa.gz") == sequences_file_path.size() - 6 )
		{

			std::vector<std::unique_ptr<SequenceFormat>> fasta_objects;
			auto fasta_parser = bioparser::createParser<bioparser::FastaParser, SequenceFormat> (sequences_file_path);

			fasta_parser -> parse_objects (fasta_objects, -1);

			statistics (fasta_objects);
		}
		else if (
			sequences_file_path.rfind(".fastq") == sequences_file_path.size() - 6 ||
                        sequences_file_path.rfind(".fq") == sequences_file_path.size() - 3  ||
                        sequences_file_path.rfind(".fastq.gz") == sequences_file_path.size() - 9 ||
                        sequences_file_path.rfind(".fq.gz") == sequences_file_path.size() - 6 )

		{
			std::vector<std::unique_ptr<SequenceFormat>> fastq_objects;
                        auto fastq_parser = bioparser::createParser<bioparser::FastqParser, SequenceFormat> (sequences_file_path);

                        fastq_parser -> parse_objects (fastq_objects, -1);

			statistics (fastq_objects);

		}

		else
			std::cout << "Unsuported format. File containing sequences needs to have one of the following extensions:\n"
				     "-> .fasta             -> .fastq                                                            \n"
                                     "-> .fa                -> .fq                                                               \n"
                                     "-> .fasta.gz          -> .fastq.gz                                                         \n"
                                     "-> .fa.gz             -> .fq.gz                                                          \n\n";
	}

	else
		std::cerr << "\nInvalid number of arguments.\n\n";

	return 0;
}
