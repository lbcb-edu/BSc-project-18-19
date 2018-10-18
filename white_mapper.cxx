#include<iostream>
#include<memory>
#include<vector>
#include "/home/krom/mapper/bioparser/include/bioparser/bioparser.hpp"
#define VERSION 0.1.0
using namespace std;

class SequenceFormat
{
	 public:
                const char* name;
                const char* sequence;
                uint32_t name_length;
                uint32_t sequence_length;
		friend void statistics (std::vector<std::unique_ptr<SequenceFormat>>);
		friend class FastaFormat;
		friend class FastqFormat;

	private:
		SequenceFormat (
                        const char* name, uint32_t name_length,
                        const char* sequence, uint32_t sequence_length)
                        {
                                this -> name = name;
                                this -> name_length = name_length;
                                this -> sequence = sequence;
                                this -> sequence_length = sequence_length;
			}
};

class FastaFormat: public SequenceFormat {
	public:
		friend class bioparser::FastaParser<FastaFormat>;

	private:
		FastaFormat (
			const char* name, uint32_t name_length,
                        const char* sequence, uint32_t sequence_length
			) : SequenceFormat (name, name_length, sequence, sequence_length)
			{
			}
};

class FastqFormat: public SequenceFormat {
	public:
        	const char* quality;
		uint32_t quality_length;
		friend class bioparser::FastqParser<FastqFormat>;

	private:
		FastqFormat (
			const char* name, uint32_t name_length,
                        const char* sequence, uint32_t sequence_length,
			const char* quality, uint32_t quality_length
			) : SequenceFormat (name, name_length, sequence, sequence_length)
			{
				this -> quality = quality;
				this -> quality_length = quality_length;
			}
};

void statistics (vector<unique_ptr<FastqFormat>> v)
{
        int min_length = 0;
        int max_length = 0;
        int total_length = 0;

        float avg_length = 0.0;

        cout << "\nNumber of sequences: " << v.size() <<"\n";

        for (auto &ptr : v)
        {
                if (ptr -> sequence_length < min_length)
                        min_length = ptr -> sequence_length;

                if (ptr -> sequence_length > max_length)
                        max_length = ptr -> sequence_length;

                total_length += ptr -> sequence_length;
        }

        cout << "Minimum sequence length: " << min_length << "\n";
        cout << "Maximum sequence length: " << max_length << "\n";
        cout << "Average sequence length: " << avg_length / v.size() << "\n\n";

}


int main (int argc, char* argv[])
{
	if (argc == 2)
	{
		string option (argv[1]);

		if (option.compare("-v") == 0)
			cout << "\nv0.1.0\n\n";

		else if (option.compare("-h") == 0)
			cout << "\nThis program is a mapper that requires two files. First one containing a set\n"
				"of fragments in FASTA or FASTQ format and the second one containing corresponding\n"
				"reference genome in FASTA format. The mapper parses the files and stores them in\n"
				"memory while also providing statistics of given files:\n"
				"-> number of sequences\n"
				"-> average length\n"
				"-> minimal and maximal length.\n\n"
				"There are two options:\n"
				"-> \"-h\" for help\n"
				"-> \"-v\" for displaying the current version.\n\n";
		else
			cerr << "\nInvalid option \"" << argv[1] << "\".\n\n";
	}

	else if (argc == 3)
	{

		string sequence_file_path (argv[1]);
		string genome_file_path (argv[2]);

		if (
			sequence_file_path.rfind(".fasta") == sequence_file_path.size() - 6 ||
			sequence_file_path.rfind(".fa") == sequence_file_path.size() - 3 ||
			sequence_file_path.rfind(".fasta.gz") == sequence_file_path.size() - 9 ||
			sequence_file_path.rfind(".fa.gz") == sequence_file_path.size() - 6 )
		{

			vector<unique_ptr<FastaFormat>> fasta_objects;
			auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FastaFormat> (sequence_file_path);

			fasta_parser -> parse_objects (fasta_objects, -1);

			int min_length = fasta_objects.at(1) -> sequence_length;
			int max_length = fasta_objects.at(1) -> sequence_length;
			int total_length = 0;

			float avg_length = 0.0;

			cout << "\nNumber of sequences: " << fasta_objects.size() <<"\n";

			for (auto &ptr : fasta_objects)
        		{
                		if (ptr -> sequence_length < min_length)
                        		min_length = ptr -> sequence_length;

                		if (ptr -> sequence_length > max_length)
                        		max_length = ptr -> sequence_length;

                		total_length += ptr -> sequence_length;
        		}

        		cout << "Minimum sequence length: " << min_length << "\n";
        		cout << "Maximum sequence length: " << max_length << "\n";
        		cout << "Average sequence length: " << (float) total_length / fasta_objects.size() << "\n\n";

		}

		else if (
			sequence_file_path.rfind(".fastq") == sequence_file_path.size() - 6 ||
                        sequence_file_path.rfind(".fq") == sequence_file_path.size() - 3  ||
                        sequence_file_path.rfind(".fastq.gz") == sequence_file_path.size() - 9 ||
                        sequence_file_path.rfind(".fq.gz") == sequence_file_path.size() - 6 )

		{
			vector<unique_ptr<FastqFormat>> fastq_objects;
                        auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FastqFormat> (sequence_file_path);

                        fastq_parser -> parse_objects (fastq_objects, -1);

                        int min_length = 0;
        		int max_length = 0;
        		int total_length = 0;

        		float avg_length = 0.0;

        		cout << "\nNumber of sequences: " << fastq_objects.size() <<"\n";

        		for (auto &ptr : fastq_objects)
        		{
                		if (ptr -> sequence_length < min_length)
                        		min_length = ptr -> sequence_length;

                		if (ptr -> sequence_length > max_length)
                        		max_length = ptr -> sequence_length;

                		total_length += ptr -> sequence_length;
        		}

        		cout << "Minimum sequence length: " << min_length << "\n";
        		cout << "Maximum sequence length: " << max_length << "\n";
        		cout << "Average sequence length: " << avg_length / fastq_objects.size() << "\n\n";
		}
	}

	else
		cerr << "\nInvalid number of arguments.\n\n";

	return 0;
}
