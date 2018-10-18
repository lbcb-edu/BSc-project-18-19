#include <iostream>
#include <vector>
#include <algorithm>
#include <bioparser/bioparser.hpp>

using namespace std;

class Fasta {
    const char* name;
    const char* sequence;
    uint32_t name_length;
    uint32_t sequence_length;
    
public:
     Fasta(const char* name, uint32_t name_length,
        const char* sequence, uint32_t sequence_length) {
            this -> name = name;
            this -> name_length = name_length;
            this -> sequence = sequence;
            this -> sequence_length = sequence_length; 
        }
};

class Fastq {
        const char* name;
        const char* sequence;
        const char* quality;
        uint32_t name_length;
        uint32_t sequence_length;
        uint32_t quality_length;
public:
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

std::vector<std::unique_ptr<Fasta>> parseFasta (string fastaFile){
    std::vector<std::unique_ptr<Fasta>> fasta_objects;
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Fasta>(&fastaFile);
    fasta_parser->parse_objects(fasta_objects, -1);
    return fasta_objects;
}

 std::vector<std::unique_ptr<Fastq>> parseFastq (string fastqFile){
    std::vector<std::unique_ptr<Fastq>> fastq_objects;
    auto fastq_parser = bioparser::createParser<bioparser::FastqParser, Fastq>(&fastqFile);
    uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
    while (true) {
        auto status = fastq_parser->parse_objects(fastq_objects, size_in_bytes);
        if (status == false) {
            break;
        }
    }
    return fastq_objects;
}

int numOfSeq (std::vector<std::unique_ptr<Fastq>> fastq_objects){
    return 
    }

bool isFasta(string arg){
    bool a;
    arg = "00000000" + arg;
    int len = arg.length();
    std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
//    cout << arg<< endl;
//    cout << arg.substr(len - 6)<<endl;
//    cout << len << endl;
//    a = arg.substr(len - 6).compare(".fasta")    == 0;
    a = arg.substr(len - 6).compare(".fasta")    == 0 ||
        arg.substr(len - 3).compare(".fa")   == 0 ||
        arg.substr(len - 9).compare(".fasta.gz") == 0 ||
        arg.substr(len - 6).compare(".fa.gz")    == 0;
    if (a) cout << "a" << endl;
    return a;
}

bool isFastq(string argg){
    argg = "        " + argg;
    int len = argg.length();
    cout << "usloq" << endl;
    std::transform(argg.begin(), argg.end(), argg.begin(), ::tolower);
    return  argg.substr(len - 6).compare(".fastq")    == 0 ||
        argg.substr(len - 3).compare(".fq")       == 0 ||
        argg.substr(len - 9).compare(".fastq.gz") == 0 ||
        argg.substr(len - 6).compare(".fq.gz")    == 0;
}

int main(int argc, char* argv[]) {
    vector<string> allArgs(argv, argv+argc);
    
    string arg1;
    string arg2;
    std::vector<std::unique_ptr<Fastq>> fastq_objects;
    std::vector<std::unique_ptr<Fasta>> fasta_objects1;
    std::vector<std::unique_ptr<Fasta>> fasta_objects2;
    
    switch(argc) {
        case 2: arg1 = allArgs.at(1);
                break;
                 
        case 3: arg1 = allArgs.at(1);
                arg2 = allArgs.at(2);
                break;
        
        default: cout << "Wrong number of arguments." << endl;
                 return 1;
    }
    
    if(arg1.compare("-h") == 0 || arg1.compare("--help") == 0) {
        cout << "Poruka!" << endl;
        
    } else if(arg1.compare("-v") == 0 || arg1.compare("--version") == 0) {
        cout << "v0.1.0" << endl;
        
    } else if(argc == 3 && (isFasta(arg1) || isFastq(arg1)) && isFasta(arg2)){
        if (isFasta(arg1)) {
            fasta_objects1 = parseFasta(arg1);
        } else {
            fastq_objects = parseFastq(arg1);
        }
        fasta_objects2 = parseFasta(arg2);
        
            
    } else {
        cout << "Wrong input." << endl;
        return 1;
    }
    
    return 0;
}  