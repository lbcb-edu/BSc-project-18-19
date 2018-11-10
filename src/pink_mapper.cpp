#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <getopt.h>
#include <cstring>
#include <stdlib.h>
#include <time.h>

#include <bioparser/bioparser.hpp>
#include "pink_alignment.hpp"

static struct option options[] = {
    {"help",    no_argument, 0, 'h'},
    {"version", no_argument, 0, 'v'},
    {NULL,      no_argument, 0,  0 }
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

std::vector<std::unique_ptr<Fast>> parseFasta (std::string fastaFile) {
    std::vector<std::unique_ptr<Fast>> fasta_objects;
        
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Fast>(fastaFile);
    fasta_parser->parse_objects(fasta_objects, -1);

    return fasta_objects;
}

std::vector<std::unique_ptr<Fast>> parseFastq (std::string fastqFile) {
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

void printStats(const std::vector<std::unique_ptr<Fast>> &fast_objects) {
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

bool correctExtension(std::string arg, char ext_flag) {
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

pink::AlignmentType alignmentType(std::string type) {
    if (type.compare("global") == 0) {
        return pink::global;
    
    } else if (type.compare("semi_global") == 0) {
        return pink::semi_global;
        
    } else if (type.compare("local") == 0) {
        return pink::local;
    
    } else {
        throw -1;
    }
}

void help() {
    printf(
        "usage: pink_mapper [options ...] <fragments> <genome> <alignment> <match> <mismatch> <cost>\n"
        "\n"
        "   <fragments>\n"
        "       input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "       containing set of fragments\n"
        "   <genome>\n"
        "       input file in FASTA format (can be compressed with gzip)\n"
        "       containing corresponding reference genome\n"
        "   <alignment>\n"
        "        input type of alignment: global, semi_global or local\n"
        "   <match>\n"
        "       input match cost\n"
        "   <mismatch>\n"
        "       input mismatch cost\n"
        "   <gap>\n"
        "       input insertion/deletion cost\n"
        "\n"
        "   options:\n"
        "       -v, --version\n"
        "           prints the version number\n"
        "       -h, --help\n"
        "           prints the usage\n");
}

void version() {
    printf("v0.1.0\n");
}

int main(int argc, char* argv[]) {
    std::string helpMessage = "Wrong input. Use \"-h\" or \"--help\" for help.";
    srand (time(NULL));
    
    if (argc - optind == 1) {
        char optchr;
        
        while((optchr = getopt_long(argc, argv, "hv", options, NULL)) != -1) {
            switch(optchr) {
                case 'h': help();
                          return 0;
                case 'v': version();
                          return 0;
                default:  std::cerr << helpMessage << std::endl;
                          return 1;
            }
        }
    }

    if (argc - optind != 2 && argc - optind != 6) {
        std::cout << helpMessage << std::endl;
        return 1;
    }

    if ((correctExtension(argv[1], 'a') || correctExtension(argv[1], 'q')) && correctExtension(argv[2], 'a')) {
        std::vector<std::unique_ptr<Fast>> fast_objects1;
        std::vector<std::unique_ptr<Fast>> fast_objects2;
        
        if (correctExtension(argv[1], 'a')) {
            fast_objects1 = parseFasta(argv[1]);
        } else {
            fast_objects1 = parseFastq(argv[1]);
        }
        
        fast_objects2 = parseFasta(argv[2]);
       
        if(argc - optind == 2) {
            std::cerr << "~FIRST FILE~" << std::endl;
            printStats(fast_objects1);
            std::cerr << "\n" << "~SECOND FILE~" << std::endl;
            printStats(fast_objects2);
            return 0;
        }
        
        int query  = rand() % fast_objects1.size();
        int target = rand() % fast_objects1.size();
        
//        const char* q = "ATCCGAT";
//        const char* t = "TGCATAT";
//        unsigned q_len = 7;
//        unsigned t_len = 7;
//
        const char* q  = (fast_objects1[query]  -> sequence).c_str();
        const char* t  = (fast_objects1[target] -> sequence).c_str();
        unsigned q_len = (fast_objects1[query]  -> sequence).length();
        unsigned t_len = (fast_objects1[target] -> sequence).length();

        pink::AlignmentType type;
        int match;
        int mismatch;
        int gap;
        std::string cigar;
        unsigned int target_begin = 0;
        
        try {
            type     = alignmentType(argv[3]);
            match    = std::stoi(argv[4]);
            mismatch = std::stoi(argv[5]);
            gap      = std::stoi(argv[6]);
        
        } catch (int e) {
            std::cout << helpMessage << std::endl;
            return 1;
        }
        
        int cost = pink::pairwise_alignment(q, q_len, t, t_len, type, match, mismatch, gap);
        std::cout << "\nFinal cost: " << cost << std::endl;
        
        pink::pairwise_alignment(q, q_len, t, t_len, type, match, mismatch, gap, cigar, target_begin);
        cigar = std::string(cigar.rbegin(), cigar.rend());
        std::cout << "\nCigar: " << cigar << "\n\n";
            
    } else {
        std::cout << helpMessage << std::endl;
        return 1;
    }
    
    return 0;
}  
