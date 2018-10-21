#include <iostream>
#include <vector>
#include <algorithm>
#include <bioparser/bioparser.hpp>

class Fasta {

public:
    const char* name;
    const char* sequence;
    uint32_t name_length;
    uint32_t sequence_length;

    Fasta(
        const char* name, uint32_t name_length,
        const char* sequence, uint32_t sequence_length) {
            this -> name = name;
            this -> name_length = name_length;
            this -> sequence = sequence;
            this -> sequence_length = sequence_length; 
        }
};

class Fastq {

public:
    const char* name;
    const char* sequence;
    const char* quality;
    uint32_t name_length;
    uint32_t sequence_length;
    uint32_t quality_length;

    Fastq(
        const char* name, uint32_t name_length,
        const char* sequence, uint32_t sequence_length,
        const char* quality, uint32_t quality_length) {
            this -> name = name;
            this -> name_length = name_length;
            this -> sequence = sequence;
            this -> sequence_length = sequence_length; 
            this -> quality = quality;
            this -> quality_length = quality_length;
    }
};

std::vector<std::unique_ptr<Fasta>> parseFasta (std::string fastaFile) {
    std::vector<std::unique_ptr<Fasta>> fasta_objects;
        
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Fasta>(fastaFile);
    fasta_parser->parse_objects(fasta_objects, -1);
        
    return fasta_objects;
}

std::vector<std::unique_ptr<Fastq>> parseFastq (std::string fastqFile) {
    std::vector<std::unique_ptr<Fastq>> fastq_objects;
        
    auto fastq_parser = bioparser::createParser<bioparser::FastqParser, Fastq>(fastqFile);
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

void printStatsFasta(std::vector<std::unique_ptr<Fasta>> fasta_objects) {
    unsigned numOfSeq = fasta_objects.size();
    unsigned sum = 0;
    float average;
    unsigned min = fasta_objects[0] -> sequence_length;
    unsigned max = min;
    
    for (unsigned i=0; i < numOfSeq;  i++) {
        sum += fasta_objects[i] -> sequence_length;
                
        if (fasta_objects[i] -> sequence_length < min) {
            min = fasta_objects[i] -> sequence_length;
        }
        if (fasta_objects[i] -> sequence_length > max) {
            max = fasta_objects[i] -> sequence_length;
        }
    }
    average = sum/numOfSeq;
    
    std::cerr << "Number of sequences: " << numOfSeq << std::endl;
    std::cerr << "Average length: "      << average  << std::endl;
    std::cerr << "Minimal length: "      << min      << std::endl;
    std::cerr << "Maximal length: "      << max      << std::endl;  
}

void printStatsFastq(std::vector<std::unique_ptr<Fastq>> fastq_objects) {
    unsigned numOfSeq = fastq_objects.size();
    unsigned sum = 0;
    float average;
    unsigned min = fastq_objects[0] -> sequence_length;
    unsigned max = min;
            
    for (unsigned i=0; i < numOfSeq;  i++) {
        sum += fastq_objects[i] -> sequence_length;
                
        if (fastq_objects[i] -> sequence_length < min) {
            min = fastq_objects[i] -> sequence_length;
        }
        if (fastq_objects[i] -> sequence_length > max) {
            max = fastq_objects[i] -> sequence_length;
        }
    }
    average = sum/numOfSeq;
    
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
    
    switch(argc) {
        case 2: arg1 = allArgs.at(1);
                break;
                 
        case 3: arg1 = allArgs.at(1);
                arg2 = allArgs.at(2);
                break;
        
        default: std::cout << "Wrong number of arguments." << std::endl;
                 return 1;
    }
    
    if(arg1.compare("-h") == 0 || arg1.compare("--help") == 0) {
        help();
        
    } else if(arg1.compare("-v") == 0 || arg1.compare("--version") == 0) {
        version();
        
    } else if(argc == 3 && (isFasta(arg1) || isFastq(arg1)) && isFasta(arg2)) {
        
        std::cerr << "~FIRST FILE~" << std::endl;
        if (isFasta(arg1)) {
            printStatsFasta(parseFasta(arg1));
        } else {
            printStatsFastq(parseFastq(arg1));
        }
        
        std::cerr << "\n" << "~SECOND FILE~" << std::endl;
        printStatsFasta(parseFasta(arg2));
        
    } else {
        std::cout << "Wrong input." << std::endl;
        return 1;
    }
    
    return 0;
}  