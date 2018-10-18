#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>

#include "OrangeConfig.h"
#include "bioparser.hpp"


using namespace std;

vector<string> FASTAExtensionVector{".fasta", ".fa", ".fasta.gz", ".fa.gz"};
vector<string> FASTQExtensionVector{".fastq", ".fq", ".fastq.gz", ".fq.gz"};

class FASTAEntity {
	
	public:
		const char *name;
		uint32_t name_length;
		const char *sequence;
		uint32_t sequence_length;
		
		FASTAEntity(
			const char *name, uint32_t name_length,
			const char *sequence, uint32_t sequence_length) {

			this->name = name;
			this->name_length = name_length;
			this->sequence = sequence;
			this->sequence_length = sequence_length;
		}
};

class FASTQEntity {

	public:
		const char *name;
		uint32_t name_length;
		const char *sequence;
		uint32_t sequence_length;
		const char* quality;
		uint32_t quality_length;
		
		FASTQEntity(
			const char* name, uint32_t name_length,
			const char* sequence, uint32_t sequence_length,
			const char* quality, uint32_t quality_length) {

			this->name = name;
			this->name_length = name_length;
			this->sequence = sequence;
			this->sequence_length = sequence_length;
			this->quality = quality;
			this->quality_length = quality_length;
		}
};

bool endsWith(string const &fullString, string const &ending) {
	if(fullString.length() < ending.length()) return false;
	
	return (fullString.compare(fullString.length() - ending.length(), ending.length(), ending) == 0);
}

bool isExtensionMemberOfVector(string const &str, vector<string> const vec) {
	for(string s : vec) {
		if(endsWith(str, s)) return true;
	}

	return false;
}

void help() {
    printf(
        "\n"
        "The program accepts two files as floating arguments and outputts number of sequences, average length, minimal and maximal length.\n"
        "Supported formats are .fasta, .fa, .fastq, .fq, .fasta.gz, .fa.gz, .fastq.gz, .fq.gz\n"
        "\n"
        "usage: orange_mapper <file_name> <file_name>\n"
        "\n"
        "options: \n"
        "      -h, --help  -  prints help menu (currently displaying)\n"
        "      -v, --version  -  prints help menu (currently displaying)\n");
}

void version() {
    printf("Version: %d\n",
        orange_mapper_VERSION_MAJOR
    );
}

void printStatistics(uint32_t max, uint32_t min, uint32_t num_of_seq, uint64_t total_length, string const &filePath) {
	fprintf(stderr, "----Statistics----\n");
	fprintf(stderr, "--%s--\n", filePath.c_str());
	fprintf(stderr, "Maximum length : %d\n", max);
	fprintf(stderr, "Minimum length : %d\n", max);
	fprintf(stderr, "Number of sequences : %d\n", num_of_seq);
	fprintf(stderr, "Average length : %f\n", (double)total_length / num_of_seq);
	fprintf(stderr, "------------------\n");
}

void readFASTQFile(string const &filePath) {
//TODO
}

void readFASTAFile(string const &filePath) {
	vector<unique_ptr<FASTAEntity>> fasta_objects;
	auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAEntity>(filePath);
	
	fasta_parser->parse_objects(fasta_objects, -1);
	
	uint32_t max = 0;

	uint32_t min = -1;

	uint32_t num_of_seq = 0;
	uint64_t total_length = 0;

	for(auto const& p : fasta_objects) {
		num_of_seq++;
		total_length += p-> sequence_length;

		max = (max > p -> sequence_length) ? max : p-> sequence_length;
			
		if(min == -1) {
			min = p-> sequence_length;
		} else {
			min = (min < p -> sequence_length) ? min : p-> sequence_length;
		}
	}

	printStatistics(max, min, num_of_seq, total_length, filePath);
}


void calculateAndPrintOutStatistics(string const &firstFilePath, string const &secondFilePath, bool isFirstFileFASTA) {
	if(isFirstFileFASTA) {
		readFASTAFile(firstFilePath);
	} else {
		readFASTQFile(firstFilePath);
	}

	readFASTAFile(secondFilePath);
}

int main(int argc, char** argv) {

	struct option options[] = {
		{"version", no_argument, 0, 'v'},
		{"help", no_argument, 0, 'h'},
	};

	char optchr;

	while((optchr = getopt_long(argc, argv, "hv", options, NULL)) != -1) {
		switch(optchr) {
			case 0:
				break;
			case 'h':
				help();
				return 0;
			case 'v':
				version();
				return 0;
			default:
				fprintf(stderr, "Entered option is not valid.\n");
				fprintf(stderr, "Use \"-h\" or \"--help\" for more information.\n");
				return 1;
		}
	}

	if(argc - 1 < 2) {
		fprintf(stderr, "Program requires more than one argument!\n");
		fprintf(stderr, "Use \"-h\" or \"--help\" for more information.\n");
		return 1;
	}

	string firstFilePath = argv[1];
	string secondFilePath = argv[2];

	bool isFirstFASTA = isExtensionMemberOfVector(firstFilePath, FASTAExtensionVector);

	if(!((isFirstFASTA || isExtensionMemberOfVector(firstFilePath, FASTQExtensionVector))
		&& isExtensionMemberOfVector(secondFilePath, FASTAExtensionVector))) {
		fprintf(stderr, "One or more given arguments are not valid file format!\n");
		fprintf(stderr, "Use \"-h\" or \"--help\" for more information.\n");
		return 1;
	}

	calculateAndPrintOutStatistics(firstFilePath, secondFilePath, isFirstFASTA);

	return 0;
}
