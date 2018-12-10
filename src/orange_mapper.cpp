#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <string.h>
#include <fstream>

#include "OrangeConfig.h"
#include "bioparser/bioparser.hpp"
#include "orange_alignment.h"
#include "orange_minimizers.h"

std::vector<std::string> FASTAExtensionVector{".fasta", ".fa", ".fasta.gz", ".fa.gz"};
std::vector<std::string> FASTQExtensionVector{".fastq", ".fq", ".fastq.gz", ".fq.gz"};

struct option options[] = {
		{"version", no_argument, 0, 'v'},
		{"help", no_argument, 0, 'h'},
		{"global", no_argument, 0, 'g'},
		{"semi-global", no_argument, 0, 's'},
		{"local", no_argument, 0, 'l'},
		{"match", required_argument, 0, 0},
		{"mismatch", required_argument, 0, 0},
		{"gap", required_argument, 0, 0},
		{"kmer", required_argument, 0, 0},
		{"window_lenght", required_argument,0 ,0},
		{0, 0, 0, 0}
	};

struct alignmentStruct {
	orange::AlignmentType type = orange::AlignmentType::no_alignment;
	int match = 1;
	int mismatch = -1;
	int gap = -1;
};
typedef struct alignmentStruct alignment;

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
		std::string name;
		std::string sequence;
		
		FASTAEntity(
			const char *name, uint32_t name_length,
			const char *sequence, uint32_t sequence_length) {

			(this->name).assign(name, name_length);

			(this->sequence).assign(sequence, sequence_length);
		}
};

class FASTQEntity : public FASTAEntity {

	public:
		std::string quality;
		
		FASTQEntity(
			const char* name, uint32_t name_length,
			const char* sequence, uint32_t sequence_length,
			const char* quality, uint32_t quality_length)
			: FASTAEntity(name, name_length, sequence, sequence_length) {
			
			(this->quality).assign(quality, quality_length);

		}
};

bool endsWith(std::string const &fullString, std::string const &ending) {
	if(fullString.length() < ending.length()) return false;
	
	return (fullString.compare(fullString.length() - ending.length(), ending.length(), ending) == 0);
}

bool isExtensionMemberOfVector(std::string const &str, std::vector<std::string> const vec) {
	for(auto const& s : vec) {
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
        "	-h, --help  -  prints help menu (currently displaying)\n"
        "	-v, --version  -  prints program version\n"
		"	-g, --global - alongside statistics, prints CIGAR string gotten as result of global alignment of 2 random sequences in the first file\n"
		"	-l, --local - alongside statistics, prints CIGAR string gotten as result of local alignment of 2 random sequences in the first file\n"
		"	-s, --semi-global - alongside statistics, prints CIGAR string gotten as result of semi-global alignment of 2 random sequences in the first file\n"
		"	--match (arg) - sets match score parameter that is used in alignments(default value is 1)\n"
		"	--mismatch (arg) - sets mismatch score parameter that is used in alignments(default value is -1)\n"
		"	--gap (arg) - sets gap score parameter that is used in alignments(default value is -1)\n\n"
		"   -k (arg) - sets size of minimizer"
		"   -w (arg) - sets window lenght for minimizers");
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
void calculateStats(std::vector<std::unique_ptr<T>> const &entities, stats *fileStats, alignment alignment) {
	
	for(auto const& p : entities) {
		fileStats->num_of_seq++;
		fileStats->total_length += (p-> sequence).length();

		fileStats->max = (fileStats->max > (p-> sequence).length()) ? fileStats->max : (p-> sequence).length();
			
		if(fileStats->min == -1) {
			fileStats->min = (p-> sequence).length();
		} else {
			fileStats->min = (fileStats->min < (p-> sequence).length()) ? fileStats->min : (p-> sequence).length();
		}
	}

	if(alignment.type != orange::AlignmentType::no_alignment) {
		int firstRand = rand()%entities.size();
		int secondRand = rand()%entities.size();		

		while(secondRand == firstRand) {
			secondRand = rand()%entities.size();
		}

		auto const &query = entities[firstRand];
		auto const &target = entities[secondRand];

		std::string cigar;
		unsigned int target_begin;

		orange::pairwise_alignment(query->sequence.c_str(), query->sequence.length(), target->sequence.c_str(), target->sequence.length(), 
		alignment.type, alignment.match, alignment.mismatch, alignment.gap, cigar, target_begin);

		printf("START POSITION: %d\n", target_begin);

		printf("--CIGAR:--\n");
		printf("%s\n\n", cigar.c_str());
	}

}

void readFASTQFile(std::string const &filePath, alignment alignment) {
	std::vector<std::unique_ptr<FASTQEntity>> fastq_objects;
	auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FASTQEntity>(filePath);

	uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
	stats fileStats;
	while (true) {
		auto status = fastq_parser->parse_objects(fastq_objects, size_in_bytes);
		if (status == false) {
			break;
		}
	}
	calculateStats<FASTQEntity>(fastq_objects, &fileStats, alignment);	
	printStatistics(fileStats, filePath);
}

void readFASTAFile(std::string const &filePath, alignment alignment) {
	std::vector<std::unique_ptr<FASTAEntity>> fasta_objects;
	auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAEntity>(filePath);
	
	fasta_parser->parse_objects(fasta_objects, -1);
	
	stats fileStats;
	calculateStats<FASTAEntity>(fasta_objects, &fileStats, alignment);

	printStatistics(fileStats, filePath);
}

void calculateAndPrintOutStatistics(std::string const &firstFilePath, std::string const &secondFilePath, bool isFirstFileFASTA, alignment alignment) {
	if(isFirstFileFASTA) {
		readFASTAFile(firstFilePath, alignment);
	} else {
		readFASTQFile(firstFilePath, alignment);
	}
	alignment.type = orange::AlignmentType::no_alignment;
	readFASTAFile(secondFilePath, alignment);
}

void createCSVfile(std::string const &filePath, int k, int window_lenght) {
	std::vector<std::unique_ptr<FASTAEntity>> fasta_objects;
	auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAEntity>(filePath);
	
	fasta_parser->parse_objects(fasta_objects, -1);

	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;
	std::map<unsigned int, unsigned int> csvMap;

	unsigned int end=fasta_objects.size();
	float counter = 0;
	float b, temp;
	float a=b=0.01;

	for (auto const &x : fasta_objects) {
		minimizers = orange::minimizers(x->sequence.c_str(), x->sequence.length(), k, window_lenght);
		counter++;
		for(auto const &y : minimizers) {
			csvMap[std::get<0>(y)]++;
		}
		if (counter/end>b) {
			temp = b*100;
			std::cout << "\b\b\b" << (int)temp << "%" << std::flush;
			b+=a;
		}
	}

	std::ofstream file;
    file.open ("minimizers.csv");

	for (auto const &m : csvMap) {
		file << m.first << "," << m.second << "\n";
	}

	std::cout << "\b\b\b Done!" << std::flush;

	file.close();

}

int main(int argc, char** argv) {

	srand((unsigned)time(0)); 
	alignment alignment;
	int k=5, window_lenght=15;

	char optchr;

	int option_index = 0;
	while((optchr = getopt_long(argc, argv, "hvgslkw", options, &option_index)) != -1) {
		switch(optchr) {
			case 0:
				if(options[option_index].flag != 0)
					break;

				if(strcmp(options[option_index].name, "match") == 0)
					alignment.match = atoi(optarg);
				else if(strcmp(options[option_index].name, "mismatch") == 0)
					alignment.mismatch = atoi(optarg);
				else if(strcmp(options[option_index].name, "gap") == 0)
					alignment.gap = atoi(optarg);

				break;
			case 'h':
				help();
				break;
			case 'v':
				version();
				break;
			case 'g':
				alignment.type = orange::AlignmentType::global;
				break;
			case 's':
				alignment.type = orange::AlignmentType::semi_global;
				break;
			case 'l':
				alignment.type = orange::AlignmentType::local;
				break;
			case 'k':
				k = atoi(optarg);
				break;
			case 'w':
				window_lenght = atoi(optarg);
				break;
			default:
				fprintf(stderr, "Entered option is not valid.\n");
				fprintf(stderr, "Use \"-h\" or \"--help\" for more information.\n");
				return 1;
		}
	}

	if(argc - optind < 2) {
		fprintf(stderr, "Program requires more than one argument!\n");
		fprintf(stderr, "Use \"-h\" or \"--help\" for more information.\n");
		return 1;
	}

	std::string firstFilePath = argv[optind];
	std::string secondFilePath = argv[optind + 1];

	bool isFirstFASTA = isExtensionMemberOfVector(firstFilePath, FASTAExtensionVector);

	if(!((isFirstFASTA || isExtensionMemberOfVector(firstFilePath, FASTQExtensionVector))
		&& isExtensionMemberOfVector(secondFilePath, FASTAExtensionVector))) {
		fprintf(stderr, "One or more given arguments are not valid file format!\n");
		fprintf(stderr, "Use \"-h\" or \"--help\" for more information.\n");
		return 1;
	}

	calculateAndPrintOutStatistics(firstFilePath, secondFilePath, isFirstFASTA, alignment);

	createCSVfile(firstFilePath, k, window_lenght);

	return 0;
}
