#include <iostream>
#include <vector>
#include <algorithm>
#include <bioparser/bioparser.hpp>
#include <getopt.h>

class Fast {

public:
    std::string name;
    std::string sequence;
    std::string quality;
    
    Fast(
        std::string name, uint32_t name_length,
        std::string sequence, uint32_t sequence_length) {
            this -> name = name;
            this -> sequence = sequence;
        }
        
    Fast(
        std::string name, uint32_t name_length,
        std::string sequence, uint32_t sequence_length,
        std::string quality, uint32_t quality_length) {
            this -> name = name;
            this -> sequence = sequence;
            this -> quality = quality;
    }
};

static struct option options[] = {
    {"help",    no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
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

bool isFasta(std::string arg) {
    arg = "        " + arg;
    int len = arg.length();
    std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
    
    return arg.substr(len - 6).compare(".fasta")    == 0 ||
           arg.substr(len - 3).compare(".fa")       == 0 ||
           arg.substr(len - 9).compare(".fasta.gz") == 0 ||
           arg.substr(len - 6).compare(".fa.gz")    == 0;
}

bool isFastq(std::string arg) {
    arg = "        " + arg;
    int len = arg.length();
    std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
        
    return arg.substr(len - 6).compare(".fastq")    == 0 ||
           arg.substr(len - 3).compare(".fq")       == 0 ||
           arg.substr(len - 9).compare(".fastq.gz") == 0 ||
           arg.substr(len - 6).compare(".fq.gz")    == 0;
}

void printStats(std::vector<std::unique_ptr<Fast>> fast_objects) {
    unsigned numOfSeq = fast_objects.size();
    unsigned sum = 0;
    float average;
    unsigned min = fast_objects[0] -> sequence.size();
    unsigned max = min;
    
    for (unsigned i=0; i < numOfSeq;  i++) {
        sum += (fast_objects[i] -> sequence).size();
                
        if ((fast_objects[i] -> sequence).size() < min) {
            min = (fast_objects[i] -> sequence).size();
        }
        if ((fast_objects[i] -> sequence.size()) > max) {
            max = (fast_objects[i] -> sequence.size());
        }
    }
    average = (float)sum / numOfSeq;
    
    std::cerr << "Number of sequences: " << numOfSeq << std::endl;
    std::cerr << "Average length: "      << average  << std::endl;
    std::cerr << "Minimal length: "      << min      << std::endl;
    std::cerr << "Maximal length: "      << max      << std::endl;  
}

void help() {
    printf(
        "usage: lbcb-mapper [options ...] <fragments> <genome>\n"
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
        "           prints the usage\n");
}

void version() {
    printf("v0.1.0\n");
}

int main(int argc, char* argv[]) {
    std::vector<std::string> allArgs(argv, argv+argc);
    
    std::string arg1;
    std::string arg2;
    
    if (argc == 2)  {
        char optchr;
        
        while((optchr = getopt_long(argc, argv, "hv", options, NULL)) != -1) {
            switch(optchr) {
                case 0: break;
                case 'h': help();
                          break;
                case 'v': version();
                          break;
                default: std:: cerr << "Wrong input. Use \"-h\" or \"--help\"" << std::endl;
                         return 1;
            }
        }
        
        if (argc  == optind) {
            return 0;
        } else if(argc - optind < 2) {
            std::cerr << "Wrong input. Use \"-h\" or \"--help\"." << std::endl;
            return 1;
        }
        
    } else if (argc == 3) {
        arg1 = allArgs.at(1);
        arg2 = allArgs.at(2);
         
        if ((isFasta(arg1) || isFastq(arg1)) && isFasta(arg2)) {
            
            std::cerr << "~FIRST FILE~" << std::endl;
            if (isFasta(arg1)) {
                printStats(parseFasta(arg1));
            } else {
                printStats(parseFastq(arg1));
            }
            
            std::cerr << "\n" << "~SECOND FILE~" << std::endl;
            printStats(parseFasta(arg2));
        }
        
    } else {
        std::cout << "Wrong input." << std::endl;
        return 1;
    }
    
    return 0;
}  
