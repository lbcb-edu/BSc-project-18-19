#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <string.h>
#include <fstream>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include <iterator>
#include <memory>

#include "OrangeConfig.h"
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "orange_alignment.h"
#include "orange_minimizers.h"

#define DEFAULT_BAND_OF_WIDTH 500

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
		{"kmer", required_argument, 0, 'k'},
		{"window_lenght", required_argument, 0 ,'w'},
		{"top_minimizers", required_argument, 0 ,'f'},
		{"cigar", no_argument, 0, 'c'},
		{"threads", required_argument, 0, 't'},
		{0, 0, 0, 0}
	};

struct alignmentStruct {
	orange::AlignmentType type = orange::AlignmentType::no_alignment;
	int match = 1;
	int mismatch = -1;
	int gap = -1;
};
typedef struct alignmentStruct alignment;

auto cmp_def = [](std::pair<unsigned int,unsigned int> const & a, std::pair<unsigned int,unsigned int> const & b) { 
 		return a.second != b.second?  a.second > b.second : a.first > b.first;
	};

auto custom_cmp_1 = [](std::tuple<short, long int, long int> const & a, std::tuple<short, long int, long int> const & b) { 
			if(std::get<0>(a) == 1) return std::get<1>(a) - std::get<2>(a) < std::get<1>(b) - std::get<2>(b);
			return std::get<1>(a) + std::get<2>(a) < std::get<1>(b) + std::get<2>(b);
	};// ???

auto custom_cmp = [](std::tuple<short, long int, unsigned int> const & a, std::tuple<short, long int, unsigned int> const & b) { 
				//if(std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);		
				if(std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
				//if(std::get<2>(a) != std::get<2>(b))
				return std::get<1>(a) < std::get<1>(b);
				//return std::get<3>(a) < std::get<3>(b);
			};// ???

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

class FASTAQEntity {
	
	public:
		std::string name;
		std::string sequence;
		std::string quality;
		
		FASTAQEntity(
			const char *name, uint32_t name_length,
			const char *sequence, uint32_t sequence_length) {

			(this->name).assign(name, name_length);

			(this->sequence).assign(sequence, sequence_length);
		}

		FASTAQEntity(
			const char* name, uint32_t name_length,
			const char* sequence, uint32_t sequence_length,
			const char* quality, uint32_t quality_length) : FASTAQEntity(name, name_length, sequence, sequence_length) {
			
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
		"	-k (arg) - sets size of minimizer"
		"	-w (arg) - sets window length for minimizers"
		"	-f (arg) - sets number of ignored minimizers in promille");
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

void calculateStats(std::vector<std::unique_ptr<FASTAQEntity>> const &entities, stats *fileStats, alignment alignment) {
	
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

std::vector<std::unique_ptr<FASTAQEntity>> readFASTQFile(std::string const &filePath, alignment alignment, bool calc_statistics) {
	std::vector<std::unique_ptr<FASTAQEntity>> fastq_objects;
	auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FASTAQEntity>(filePath);

	uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
	while (true) {
		auto status = fastq_parser->parse_objects(fastq_objects, size_in_bytes);
		if (status == false) {
			break;
		}
	}
	
	if(calc_statistics) {
		stats fileStats;
		calculateStats(fastq_objects, &fileStats, alignment);	
		printStatistics(fileStats, filePath);
	}

	return fastq_objects;
}

std::vector<std::unique_ptr<FASTAQEntity>> readFASTAFile(std::string const &filePath, alignment alignment, bool calc_statistics) {
	std::vector<std::unique_ptr<FASTAQEntity>> fasta_objects;
	auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAQEntity>(filePath);
	
	fasta_parser->parse_objects(fasta_objects, -1);
	
	if(calc_statistics) {
		stats fileStats;
		calculateStats(fasta_objects, &fileStats, alignment);
		printStatistics(fileStats, filePath);
	}

	return fasta_objects;
}

void calculateAndPrintOutStatistics(std::string const &firstFilePath, std::string const &secondFilePath, bool isFirstFileFASTA, alignment alignment) {
	if(isFirstFileFASTA) {
		readFASTAFile(firstFilePath, alignment, true);
	} else {
		readFASTQFile(firstFilePath, alignment, true);
	}
	alignment.type = orange::AlignmentType::no_alignment;
	readFASTAFile(secondFilePath, alignment, true);
}

void findMinimizers(std::string const &filePath, int k, int window_lenght, double f, bool isFirstFASTA) {
	alignment al;
	std::vector<std::unique_ptr<FASTAQEntity>> fastaq_objects = 
		isFirstFASTA ? readFASTAFile(filePath, al, false) : readFASTQFile(filePath, al, false);

	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;
	std::unordered_map<unsigned int, unsigned int> csvMap;

	unsigned int end=fastaq_objects.size();
	unsigned int counter = 0;

	std::cout << "Looking for minimizers, please wait\n";

	for (auto const &x : fastaq_objects) {
		minimizers = orange::minimizers(x->sequence.c_str(), x->sequence.length(), k, window_lenght);

		++counter;

		std::cout << "\r" << "Progress: " << counter << "/" << end << std::flush;

		for(auto const &y : minimizers) {
			++csvMap[std::get<0>(y)];
		}
	}

	printf("\nNumber of minimizers: %lu\n", csvMap.size());

	std::vector<std::pair<unsigned int, unsigned int>> vec(csvMap.begin(), csvMap.end());

	std::sort(vec.begin(), vec.end(), cmp_def);

	int pos = f*csvMap.size();

	printf("%d", vec[pos].second);
	
	std::cout << "\nDone!\n";
}

std::unordered_map<unsigned int, std::vector<std::tuple<unsigned int, bool>>> constructMinimizerIndex(double f, int k, int window_lenght, std::unique_ptr<FASTAQEntity> const &reference_gen) {
	printf("Looking for reference minimizers...\n");

	std::vector<std::tuple<unsigned int, unsigned int, bool>> ref_minimizers = orange::minimizers(reference_gen->sequence.c_str(), reference_gen->sequence.length(), k, window_lenght);
	printf("Reference minimizers found!\n");

	printf("Creating reference index...\n");
	std::unordered_map<unsigned int, std::vector<std::tuple<unsigned int, bool>>> ref_index;

	for(auto const &t : ref_minimizers) {
		ref_index[std::get<0>(t)].emplace_back(std::get<1>(t), std::get<2>(t));
	}
	printf("Reference index created!\n");

	printf("Removing too frequent minimizers from reference index...\n");
	auto cmp = [](std::pair<unsigned int, std::vector<std::tuple<unsigned int, bool>>> const & a, std::pair<unsigned int, std::vector<std::tuple<unsigned int, bool>>> const & b) { 
 		return a.second.size() > b.second.size();
	};
	std::vector<std::pair<unsigned int, std::vector<std::tuple<unsigned int, bool>>>> temp_vec(ref_index.begin(), ref_index.end());

	std::sort(temp_vec.begin(), temp_vec.end(), cmp);

	unsigned int border = f * temp_vec.size();
	for(unsigned int i = 0; i < border; i++) {
		ref_index.erase(temp_vec[i].first);
	}
	printf("Too frequent reference minimizers removed!\n");

	return ref_index;
}

void LISAlgorithm(std::vector<std::tuple<short, long int, long int>> const &vec, unsigned int start, unsigned int end, unsigned int &query_start, unsigned int &query_end, unsigned int &ref_start, unsigned int &ref_end, unsigned int &max_len) {
	int pArr[end - start + 1];
	int mArr[end - start + 2];

	int l = 0;

	for(int i = 0; i < end - start + 1; ++i) {
		int lo = 1;
		int hi = l;

		while(lo <= hi) {
			int mid = (lo + hi)/2;
			if(std::get<2>(vec[start + mArr[mid]]) < std::get<2>(vec[start + i]))
				lo = mid + 1;
			else
				hi = mid - 1;

		}

		int newL = lo;

		pArr[i] = mArr[newL - 1];
		mArr[newL] = i;

		if(newL > l) l = newL;
	}

	if(max_len > l) return;

	max_len = l;

	int k = mArr[l];

	for(int i = l - 1; i >= 0; --i) {
		if(i == l - 1) {
			query_end = std::get<0>(vec[start + k]) == 0 ? std::get<1>(vec[start + k]) + std::get<2>(vec[start + k]) : std::get<1>(vec[start + k]) - std::get<2>(vec[start + k]);
			ref_end = std::get<2>(vec[start + k]);
		}

		if(i == 0) {
			query_start = std::get<0>(vec[start + k]) == 0 ? std::get<1>(vec[start + k]) + std::get<2>(vec[start + k]) : std::get<1>(vec[start + k]) - std::get<2>(vec[start + k]);
			ref_start = std::get<2>(vec[start + k]);
		}
	
		k = pArr[k];
	}
}

void findLongestLinearChain(std::vector<std::tuple<short, long int, long int>> &vec, unsigned int start, unsigned int end, unsigned int &query_start, unsigned int &query_end, unsigned int &ref_start, unsigned int &ref_end, unsigned int &max_len) {
	std::sort(vec.begin() + start, vec.begin() + end + 1, custom_cmp_1);

	//for(int i  = 0; i < vec.size(); i++) {
	//	printf("++ %d, %d, %d ++ \n", std::get<0>(vec[i]), std::get<1>(vec[i]), std::get<2>(vec[i]));
	//}
	//printf("\n\n");


	LISAlgorithm(vec, start, end, query_start, query_end, ref_start, ref_end, max_len);
}

void mapThread (std::vector<std::unique_ptr<FASTAQEntity>> const &fastaq_objects, unsigned int k, unsigned int window_lenght, std::unordered_map<unsigned int, std::vector<std::tuple<unsigned int, bool>>> const &ref_index, std::unique_ptr<FASTAQEntity> const &y, int gap, int mismatch, int match, int threads, int j) {

	unsigned int target_begin;
	for(unsigned int i = j; i < fastaq_objects.size(); i += threads) {
		std::vector<std::tuple<unsigned int, unsigned int, bool>> fragment_minimizers = orange::minimizers(fastaq_objects[i]->sequence.c_str(), fastaq_objects[i]->sequence.length(), k,window_lenght);
		std::vector<std::tuple<short, long int, long int>> vec;
		for(auto const &t : fragment_minimizers) {
			if(ref_index.find(std::get<0>(t)) == ref_index.end()) continue;
			for(auto const &rt : ref_index.at(std::get<0>(t))) {
				if(std::get<2>(t) == std::get<1>(rt)) {
					vec.emplace_back(0, std::get<1>(t) - std::get<0>(rt), std::get<0>(rt));
				} else {
					vec.emplace_back(1, std::get<1>(t) + std::get<0>(rt), std::get<0>(rt));
				}
			}
		}

		if(vec.size() == 0) continue;//no matches found

		std::sort(vec.begin(), vec.end(), custom_cmp);

		//for(int i  = 0; i < vec.size(); i++) {
		//	printf("++ %d, %d, %d ++ \n", std::get<0>(vec[i]), std::get<1>(vec[i]), std::get<2>(vec[i]));
		//}
		//printf("\n\n");

		unsigned int b = 0;

		unsigned int query_start;
		unsigned int query_end ;
		unsigned int ref_start;
		unsigned int ref_end;

		unsigned int max_len = 0;
		for(unsigned int i = 0; i < vec.size(); ++i) {
			if(i == vec.size() - 1 ||
				std::get<0>(vec[i + 1]) != std::get<0>(vec[i]) ||
				std::abs(std::get<1>(vec[i + 1]) - std::get<1>(vec[i])) >= DEFAULT_BAND_OF_WIDTH) {
				//printf("b = %d, i = %d\n", b, i);
				findLongestLinearChain(vec, b, i, query_start, query_end, ref_start, ref_end, max_len);
				b = i + 1;
				//printf("maxlen = %d\n", max_len);
			}
		}

		std::string cigar;
		std::string sub;
		//printf("%d %d %d %d %d\n", query_start, query_end, ref_start, ref_end, fastaq_objects[i]->sequence.length());
		return;
		orange::pairwise_alignment(fastaq_objects[i]->sequence.c_str() + query_start, query_end + k - query_start, y->sequence.c_str() + ref_start, ref_end + k - ref_start, orange::AlignmentType::global, match, mismatch, gap, cigar, target_begin);
		continue;
		unsigned int len = cigar.length();
		unsigned int count=0, e=0, sum=0;
		for(int i = 0; i < len; ++i) {
			if(cigar[i]=='I' || cigar[i]=='D') {
				sub = cigar.substr(e, i-e);
				e= i+1;
				sum += std::stoi(sub);
			}
			else {
				if(cigar[i]=='X' || cigar[i]=='=') {
					sub = cigar.substr(e, i-e);
					e = i+1;
					count += std::stoi(sub);
					sum += std::stoi(sub);
				}
			}
		}

		std::string paf;
		paf += fastaq_objects[i]->name + "\n";
		paf += std::to_string(fastaq_objects[i]->sequence.length()) + "\n";
		paf += std::to_string(query_start) + "\n";
		paf += std::to_string(query_end+k) + "\n";
		paf += "+\n";
		paf += y->name + "\n";
		paf += std::to_string(y->sequence.length()) + "\n";
		paf += std::to_string(ref_start) + "\n";
		paf += std::to_string(ref_end+k) + "\n";
		paf += std::to_string(count) + "\n";
		paf += std::to_string(sum) + "\n";
		paf += "255\n";
		paf += "cg:Z:" + cigar + "\n"; 
		printf("%s\n", paf.c_str());
	}
}

void constructAndPrintPAF(std::string const &firstFilePath, std::string const &secondFilePath, bool isFirstFASTA, int k, int window_lenght, double f, bool c, 
							int match, int mismatch, int gap, int threads) {
	alignment al;
	std::vector<std::unique_ptr<FASTAQEntity>> fastaq_objects = 
		isFirstFASTA ? readFASTAFile(firstFilePath, al, false) : readFASTQFile(firstFilePath, al, false);

	std::vector<std::unique_ptr<FASTAQEntity>> reference_gen_vec = readFASTAFile(secondFilePath, al, false);

	for(auto &y : reference_gen_vec) {
		std::unordered_map<unsigned int, std::vector<std::tuple<unsigned int, bool>>> ref_index = constructMinimizerIndex(f, k, window_lenght, y);

		std::shared_ptr<thread_pool::ThreadPool> thread_pool1 = thread_pool::createThreadPool(threads);
		std::vector<std::future<void>> thread_futures1;

		for(unsigned i = 0; i < threads; i++) {
			thread_futures1.emplace_back(thread_pool1->submit_task(mapThread, std::ref(fastaq_objects), k, window_lenght, std::ref(ref_index), std::ref(y), gap, mismatch, match, threads, i));
		}

		for (auto &it : thread_futures1) {
			it.wait();
		}

		//mapThread(fastaq_objects, k , window_lenght, ref_index, y, gap, mismatch, match);

	}

}

int main(int argc, char** argv) {

	srand((unsigned)time(0)); 
	alignment alignment;
	int k=15, window_lenght= 5, threads = 4;
	double f = 0.001;

	bool includeCIGARInPAF = false;

	char optchr;

	int option_index = 0;
	while((optchr = getopt_long(argc, argv, "hvcgslk:w:f:t:", options, &option_index)) != -1) {
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
			case 'f':
				f = std::atof(optarg);
				break;
			case 'c':
				includeCIGARInPAF = true;
				break;
			case 't':
				threads = atoi(optarg);
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

	//findMinimizers(firstFilePath, k, window_lenght, f, isFirstFASTA);

	constructAndPrintPAF(firstFilePath, secondFilePath, isFirstFASTA, k, window_lenght, f, includeCIGARInPAF, alignment.match, alignment.mismatch, alignment.gap, threads);

	return 0;
}
