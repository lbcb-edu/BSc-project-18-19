#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <getopt.h>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <unordered_map>
#include <fstream>

#include <bioparser/bioparser.hpp>
#include "pink_alignment.hpp"
#include "pink_minimizers.hpp"

static struct option options[] = {
    {"help",          no_argument,       0, 'h'},
    {"version",       no_argument,       0, 'v'},
    {"global",        no_argument,       0, 'G'},
    {"semi_global",   no_argument,       0, 'S'},
    {"local",         no_argument,       0, 'L'},
    {"match",         required_argument, 0, 'm'},
    {"mismatch",      required_argument, 0, 's'},
    {"gap",           required_argument, 0, 'g'},
	{"k",             required_argument, 0, 'k'},
	{"window_length", required_argument, 0, 'w'},
    {NULL,            no_argument,       0,  0 }
};

class Fast {
public:
    std::string name;
    std::string sequence;
    std::string quality;
    
    Fast(
        const char* name, uint32_t name_length,
        const char* sequence, uint32_t sequence_length) : 
            name {std::string (name, name_length)},
            sequence {std::string (sequence, sequence_length)} 
        {}
        
    Fast(
        const char* name, uint32_t name_length,
        const char* sequence, uint32_t sequence_length,
        const char* quality, uint32_t quality_length) :
            name {std::string (name, name_length)},
            sequence {std::string (sequence, sequence_length)}, 
            quality {std::string (quality, quality_length)}
        {}
};

std::vector<std::unique_ptr<Fast>> parse_fasta (std::string fastaFile) {
    std::vector<std::unique_ptr<Fast>> fasta_objects;
        
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Fast>(fastaFile);
    fasta_parser->parse_objects(fasta_objects, -1);

    return fasta_objects;
}

std::vector<std::unique_ptr<Fast>> parse_fastq (std::string fastqFile) {
    std::vector<std::unique_ptr<Fast>> fastq_objects;
        
    auto fastq_parser = bioparser::createParser<bioparser::FastqParser, Fast>(fastqFile);
    uint64_t size_in_bytes = 500 * 1024 * 1024;  // 500 MB
        
    while (true) {
        auto status = fastq_parser->parse_objects(fastq_objects, size_in_bytes);
        if (status == false) {
            break;
        }
    }

    return fastq_objects;
}

void print_stats(const std::vector<std::unique_ptr<Fast>> &fast_objects) {
    unsigned numOfSeq = fast_objects.size();
    unsigned sum = 0;
    float average;
    unsigned min = (fast_objects[0] -> sequence).length();
    unsigned max = min;
    
    for (unsigned i=0; i < numOfSeq; i++) {
        sum += (fast_objects[i] -> sequence).length();
                
        if ((fast_objects[i] -> sequence).length() < min) {
            min = (fast_objects[i] -> sequence).length();
        }
        if ((fast_objects[i] -> sequence).length() > max) {
            max = (fast_objects[i] -> sequence).length();
        }
    }
    average = (float)sum / numOfSeq;
    
    std::cerr << "Number of sequences: " << numOfSeq << std::endl;
    std::cerr << "Average length: "      << average  << std::endl;
    std::cerr << "Minimal length: "      << min      << std::endl;
    std::cerr << "Maximal length: "      << max      << std::endl;  
}

bool check_extension(std::string arg, char ext_flag) {
    std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
    
    std::vector<std::string> ext;
    if (ext_flag == 'a') {
        ext = {".fa", ".fasta", ".fa.gz", ".fasta.gz"};
    } else {
        ext = {".fq", ".fastq", ".fq.gz", ".fastq.gz"};
    }
                                                                                        
    for (unsigned i=0; i < ext.size(); i++) {
        if (arg.size() < ext[i].size()) {
            break;
        
        } else {
            if (arg.substr(arg.size() - ext[i].size()).compare(ext[i]) == 0) {
                return true;
            }
        }
    }
    
    return false;
}

void printError() {
	std::cerr << "Wrong input. Use \"-h\" or \"--help\" for help." << std::endl;
	exit(1);
}

int checkInput(char* optarg) {
	if(atoi(optarg) == 0 && strcmp(optarg, "0") != 0) {
		printError();
	}
	return atoi(optarg);
}

void help() {
    printf(
        "usage: pink_mapper [options ...] <fragments> <genome> \n"
        "\n"
        "   <fragments>\n"
        "       input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "       containing set of fragments\n"
        "   <genome>\n"
        "       input file in FASTA format (can be compressed with gzip)\n"
        "       containing corresponding reference genome\n"
        "\n"
        "   options:\n"
        "       -v, --version\n"
        "           prints the version number\n"
        "       -h, --help\n"
        "           prints the usage\n"
        "       -G, --global\n"
        "           set to global alignment, default alignment\n"
        "       -S, --semi_global\n"
        "           set to semi-global alignment\n"
        "       -L, --local\n"
        "           set to local alignment\n"
        "       -m, --match\n"
        "           input match cost, default: 2\n"
        "       -s, --mismatch\n"
        "           input mismatch cost, default: -1\n"
        "       -g, --gap\n"
        "           input insertion/deletion cost, default: -2\n"
		"       -k\n"
		"           input k, default: 15\n"
		"       -w\n"
		"           input window length, default: 5\n"
        );
}

void version() {
    printf("v0.1.0\n");
}

int main(int argc, char* argv[]) {
    char optchr;
    srand (time(NULL));
    
    pink::AlignmentType type = pink::global;
    int match = 2;
    int mismatch = -1;
    int gap = -2;
    std::string cigar;
    unsigned int target_begin = 0;

	unsigned int k = 15;
	unsigned int window_length = 5;
    
    while((optchr = getopt_long(argc, argv, "hvGSLm:s:g:k:w:", options, NULL)) != -1) {
        switch(optchr) {
            case 'h': 
                help();
                return 0;
            case 'v': 
                version();
                return 0;
            case 'G':
                type = pink::global;
                break;
            case 'S':
                type = pink::semi_global;
                break;
            case 'L':
                type = pink::local;
                break;
            case 'm':
				match = checkInput(optarg);
                break;
            case 's':
                mismatch = checkInput(optarg);
                break;
            case 'g':
                gap = checkInput(optarg);
                break;
			case 'k':
				k = checkInput(optarg);
				break;
			case 'w':
				window_length = checkInput(optarg);
				break;
            default:  
                printError();
        }
    }
    
    if (argc - optind != 2) {
        printError();
    }
    
    if ((check_extension(argv[optind], 'a') || check_extension(argv[optind], 'q')) && check_extension(argv[optind+1], 'a')) {
        std::vector<std::unique_ptr<Fast>> fast_objects1;
        std::vector<std::unique_ptr<Fast>> fast_objects2;
        
        if (check_extension(argv[optind], 'a')) {
            fast_objects1 = parse_fasta(argv[optind]);
        } else {
            fast_objects1 = parse_fastq(argv[optind]);
        }
        
        fast_objects2 = parse_fasta(argv[optind+1]);
       
        std::cerr << "~FIRST FILE~" << std::endl;
        print_stats(fast_objects1);
        std::cerr << "\n" << "~SECOND FILE~" << std::endl;
        print_stats(fast_objects2);
        
        
        int query  = rand() % fast_objects1.size();
        int target = rand() % fast_objects1.size();
        
        const char* q  = (fast_objects1[query]  -> sequence).c_str();
        const char* t  = (fast_objects1[target] -> sequence).c_str();
        unsigned int q_len = (fast_objects1[query]  -> sequence).length();
        unsigned int t_len = (fast_objects1[target] -> sequence).length();
        
        int cost = pink::pairwise_alignment(q, q_len, t, t_len, type, match, mismatch, gap);
        std::cout << "\nFinal cost: " << cost << std::endl;
        
        pink::pairwise_alignment(q, q_len, t, t_len, type, match, mismatch, gap, cigar, target_begin);
        cigar = std::string(cigar.rbegin(), cigar.rend());
        std::cout << "\nCigar: " << cigar << "\n\n";
        
        
		std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector;
        std::unordered_map<unsigned int, unsigned int> occurences;
        std::unordered_map<unsigned int, unsigned int>::iterator it;

		for (auto const& object: fast_objects1) {
			minimizers_vector = pink::minimizers((object -> sequence).c_str(), (object -> sequence).length(), k, window_length);

            for (auto const& minimizer: minimizers_vector) {
                it = occurences.find(std::get<0>(minimizer)); 
                
                if (it == occurences.end()) {
                    occurences.insert (std::pair<unsigned int, unsigned int>(std::get<0>(minimizer), 1));
                
                } else {
                    ++(it->second);
                }
            }
		}
        
        std::ofstream file;
        file.open ("pink_minimizers.csv");
        file << "minimizer,occurence" << "\n";
        
        for (auto o : occurences) {
            file << o.first << "," << o.second << "\n";
        }
		
        file.close();
        
        std::cout << "Work with minimizers done!" << std::endl;
        
    } else {
        printError();
    }
    
    return 0;
}  
