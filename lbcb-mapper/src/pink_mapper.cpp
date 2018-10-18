#include <iostream>
#include <vector>
#include <algorithm>
#include <bioparser/bioparser.hpp>

using namespace std;

class Fasta {
    
public:
    const char* name;
    const char* sequence;
    uint32_t name_length;
    uint32_t sequence_length;
    
     Fasta(const char* name, uint32_t name_length,
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

std::vector<std::unique_ptr<Fasta>> parseFasta (string fastaFile){
    std::vector<std::unique_ptr<Fasta>> fasta_objects;
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Fasta>(fastaFile);
    fasta_parser->parse_objects(fasta_objects, -1);
    return fasta_objects;
}

 std::vector<std::unique_ptr<Fastq>> parseFastq (string fastqFile){
    std::vector<std::unique_ptr<Fastq>> fastq_objects;
    auto fastq_parser = bioparser::createParser<bioparser::FastqParser, Fastq>(fastqFile);
    uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
    while (true) {
        auto status = fastq_parser->parse_objects(fastq_objects, size_in_bytes);
        if (status == false) {
            break;
        }
    }
    return fastq_objects;
}

bool isFasta(string arg){
    arg = "00000000" + arg;
    int len = arg.length();
    std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
    return arg.substr(len - 6).compare(".fasta")    == 0 ||
        arg.substr(len - 3).compare(".fa")   == 0 ||
        arg.substr(len - 9).compare(".fasta.gz") == 0 ||
        arg.substr(len - 6).compare(".fa.gz")    == 0;
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
    unsigned numOfSeq1;
    unsigned numOfSeq2;
    unsigned sum1 = 0;
    unsigned sum2 = 0;
    float average1 = 0;
    float average2 = 0;
    unsigned min1 = 0;
    unsigned min2 = 0;
    unsigned max1 = 0;
    unsigned max2 = 0;
    
    
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
            
            numOfSeq1 = fasta_objects1.size();
            
            max1 = fasta_objects1[0] -> sequence_length;
            for (unsigned i=0; i < fasta_objects1.size(); i++) {
                sum1 += fasta_objects1[i] -> sequence_length;
                if (fasta_objects1[i] -> sequence_length < min1) min1 = fasta_objects1[i] -> sequence_length;
                if (fasta_objects1[i] -> sequence_length > max1) max1 = fasta_objects1[i] -> sequence_length;
            }
            
        } else {
            fastq_objects = parseFastq(arg1);
            numOfSeq1 = fastq_objects.size();
            
            max1 = fastq_objects[0] -> sequence_length;
             for (unsigned i=0; i < fastq_objects.size(); i++) {
                sum1 += fastq_objects[i] -> sequence_length;
                if (fastq_objects[i] -> sequence_length < min1) min1 = fastq_objects[i] -> sequence_length;
                if (fastq_objects[i] -> sequence_length > max1) max1 = fastq_objects[i] -> sequence_length;
            }
        }
        fasta_objects2 = parseFasta(arg2);
        numOfSeq2 = fasta_objects2.size();
        
        max2 = fasta_objects2[0] -> sequence_length;
         for (unsigned i=0; i < fasta_objects2.size(); i++) {
                sum2 += fasta_objects2[i] -> sequence_length;
                if (fasta_objects2[i] -> sequence_length < min2) min2 = fasta_objects2[i] -> sequence_length;
                if (fasta_objects2[i] -> sequence_length > max2) max2 = fasta_objects2[i] -> sequence_length;
        }
            
       average1 = sum1/numOfSeq1;
       average2 = sum2/numOfSeq2;
        
        
        
            
    } else {
        cout << "Wrong input." << endl;
        return 1;
    }
    
    return 0;
}  