#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <getopt.h>
#include <bioparser/bioparser.hpp>
#include <pink_alignment/pink_alignment.hpp>

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

void printStats(const std::vector<std::unique_ptr<Fast>> &fast_objects) {
    unsigned numOfSeq = fast_objects.size();
    unsigned sum = 0;
    float average;
    unsigned min = (fast_objects[0] -> sequence).length();
    unsigned max = min;
    
    for (unsigned i=0; i < numOfSeq;  i++) {
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

void parseFasta (std::string fastaFile) {
    std::vector<std::unique_ptr<Fast>> fasta_objects;
        
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Fast>(fastaFile);
    fasta_parser->parse_objects(fasta_objects, -1);
        
    printStats(fasta_objects);
}

void parseFastq (std::string fastqFile) {
    std::vector<std::unique_ptr<Fast>> fastq_objects;
        
    auto fastq_parser = bioparser::createParser<bioparser::FastqParser, Fast>(fastqFile);
    uint64_t size_in_bytes = 500 * 1024 * 1024;  // 500 MB
        
    while (true) {
        auto status = fastq_parser->parse_objects(fastq_objects, size_in_bytes);
        if (status == false) {
            break;
        }
    }
        
    printStats(fastq_objects);
}

bool correctExtension(std::string arg, std::vector<std::string> ext) {
    std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
    
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
    std::string helpMessage = "Wrong input. Use \"-h\" or \"--help\" for help.";
    
    if (argc - optind == 1) {
        char optchr;
        
        while((optchr = getopt_long(argc, argv, "hv", options, NULL)) != -1) {
            switch(optchr) {
                case 'h': help();
                          return 0;
                case 'v': version();
                          return 0;
                default: std::cerr << helpMessage << std::endl;
                         return 1;
            }
        }
    }   
    
    if (argc - optind != 2) {
        std::cerr << helpMessage << std::endl;
        return 1;
    }
    
    std::string arg1 = argv[1];
    std::string arg2 = argv[2];
    std::vector<std::string> fastaExtensions = {".fa", ".fasta", ".fa.gz", ".fasta.gz"};
    std::vector<std::string> fastqExtensions = {".fq", ".fastq", ".fq.gz", ".fastq.gz"};

    if ((correctExtension(arg1, fastaExtensions) || correctExtension(arg1, fastaExtensions)) 
            && correctExtension(arg2, fastaExtensions)) {
            
        std::cerr << "~FIRST FILE~" << std::endl;
        if (correctExtension(arg1, fastaExtensions)) {
            parseFasta(arg1);
        } else {
            parseFastq(arg1);
        }
            
        std::cerr << "\n" << "~SECOND FILE~" << std::endl;
        parseFasta(arg2);
            
    } else {
        std::cout << helpMessage << std::endl;
        return 1;
    }
    
    return 0;
}  
