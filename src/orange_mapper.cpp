#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>

#include "OrangeConfig.h"
#include "bioparser/bioparser.hpp"

std::vector<std::string> FASTAExtensionVector{".fasta", ".fa", ".fasta.gz", ".fa.gz"};
std::vector<std::string> FASTQExtensionVector{".fastq", ".fq", ".fastq.gz", ".fq.gz"};

struct option options[] = {
		{"version", no_argument, 0, 'v'},
		{"help", no_argument, 0, 'h'},
	};

struct statsStruct {
	uint32_t max;
	uint32_t min;

	uint32_t num_of_seq;
	uint64_t total_length;

	statsStruct() {
		max = 0;
		min = -1;

		num_of_seq = 0;
		total_length = 0;
	}
};
typedef struct statsStruct stats;

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

bool endsWith(std::string const &fullString, std::string const &ending) {
	if(fullString.length() < ending.length()) return false;
	
	return (fullString.compare(fullString.length() - ending.length(), ending.length(), ending) == 0);
}

bool isExtensionMemberOfVector(std::string const &str, std::vector<std::string> const vec) {
	for(std::string s : vec) {
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
    printf("v%d.%d.0\n",
        orange_mapper_VERSION_MAJOR,
		orange_mapper_VERSION_MINOR
    );
}

void printStatistics(stats const &fileStats, std::string const &filePath) {
	fprintf(stderr, "----Statistics----\n");
	fprintf(stderr, "--%s--\n", filePath.c_str());
	fprintf(stderr, "Maximum length : %d\n", fileStats.max);
	fprintf(stderr, "Minimum length : %d\n", fileStats.min);
	fprintf(stderr, "Number of sequences : %d\n", fileStats.num_of_seq);
	fprintf(stderr, "Average length : %f\n", (double)fileStats.total_length / fileStats.num_of_seq);
	fprintf(stderr, "------------------\n");
}

template<class T>
void calculateStats(std::vector<std::unique_ptr<T>> const &entities, stats *fileStats) {
	
	for(auto const& p : entities) {
		fileStats->num_of_seq++;
		fileStats->total_length += p-> sequence_length;

		fileStats->max = (fileStats->max > p -> sequence_length) ? fileStats->max : p-> sequence_length;
			
		if(fileStats->min == -1) {
			fileStats->min = p-> sequence_length;
		} else {
			fileStats->min = (fileStats->min < p -> sequence_length) ? fileStats->min : p-> sequence_length;
		}
	}
}

void readFASTQFile(std::string const &filePath) {
	std::vector<std::unique_ptr<FASTQEntity>> fastq_objects;
	auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FASTQEntity>(filePath);

	uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
	stats fileStats;
	while (true) {
		auto status = fastq_parser->parse_objects(fastq_objects, size_in_bytes);
		calculateStats<FASTQEntity>(fastq_objects, &fileStats);	
		
		if (status == false) {
			break;
		}
	}

	printStatistics(fileStats, filePath);
}

void readFASTAFile(std::string const &filePath) {
	std::vector<std::unique_ptr<FASTAEntity>> fasta_objects;
	auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAEntity>(filePath);
	
	fasta_parser->parse_objects(fasta_objects, -1);
	
	stats fileStats;
	calculateStats<FASTAEntity>(fasta_objects, &fileStats);

	printStatistics(fileStats, filePath);
}

void calculateAndPrintOutStatistics(std::string const &firstFilePath, std::string const &secondFilePath, bool isFirstFileFASTA) {
	if(isFirstFileFASTA) {
		readFASTAFile(firstFilePath);
	} else {
		readFASTQFile(firstFilePath);
	}

	readFASTAFile(secondFilePath);
}

int main(int argc, char** argv) {

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

	std::string firstFilePath = argv[1];
	std::string secondFilePath = argv[2];

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
