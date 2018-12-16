#include <iostream>
#include <fstream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <memory>
#include <stdint.h>
#include <getopt.h>
#include <map>
#include <unordered_map>
#include "bioparser/bioparser.hpp"
#include "blue_alignment.hpp"
#include "blue_minimizers.hpp"

namespace std {
template <> struct hash<std::tuple<unsigned int, unsigned int, bool >> {
    inline size_t operator()(const std::tuple<unsigned int, unsigned int, bool > &v) const {
        std::hash<int> int_hasher;
        return int_hasher(std::get<0>(v)) ^ int_hasher(std::get<1>(v)) ^ int_hasher(std::get<2>(v));
    }
};
}


class InputFile {
    public:
        std::string name;
        std::string sequence;
        std::string quality;

     InputFile(
       const char* name_, uint32_t name_length_,
       const char* sequence_, uint32_t sequence_length_,
       const char* quality_, uint32_t quality_length_
    ) :
         name(name_, name_length_),
         sequence(sequence_, sequence_length_),
         quality(quality_, quality_length_)

     { }

     InputFile(
       const char* name_, uint32_t name_length_,
        const char* sequence_, uint32_t sequence_length_
    ) :
         name(name_, name_length_),
         sequence(sequence_, sequence_length_)
     { }
};

static struct option long_options[] = {
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    {"match", required_argument, NULL, 'm' },
    {"mismatch", required_argument, NULL, 's'},
    {"gap", required_argument, NULL, 'g'},
    {"type", required_argument, NULL, 't'},
    {"kmer_length", required_argument, NULL, 'k'},
    {"window_length", required_argument, NULL, 'w'},
    {NULL, no_argument, NULL, 0}
};

template<typename T>
void fastaq_stat(std::vector<std::unique_ptr<T>>& fq_objects) {

    uint32_t max=0, min=INT_MAX;
    uint64_t sum=0;

    for (auto& i : fq_objects) {
        uint32_t my_length = (i->sequence).length();
        if (my_length>max){
            max=my_length;
        }

        if (my_length<min) {
            min=my_length;
        }

        sum+=my_length;
    }

    std::cout << "Maximum length of a sequence is: " << max << std::endl;
    std::cout << "Minimum length of a sequence is: " << min << std::endl;
    std::cout << "Average length of a sequence is: " << sum/(double)fq_objects.size() << std::endl;
    std::cout << "Number of sequences is: "<< fq_objects.size() << '\n' << std::endl;
};

int main (int argc, char* argv[]) {

        int c;
        int match = 0;
        int mismatch = -1;
        int gap = -1;
        char* align_type = "global";

        unsigned int kmer_length = 5;
        unsigned int window_length = 15;

        while ((c = getopt_long (argc, argv, "hvm:s:g:t:k:w:", long_options, NULL)) != -1) {
            switch(c) {
                case 'h':
                    std::cout << "You've asked for help.\n\nThis program implements 3 algorithms for pairwise alignment:" << std::endl
                              << "\t- Needleman-Wunsch algorithm for global alignment" << std::endl
                              << "\t- Smith-Waterman algorithm for local alignment" << std::endl
                              << "\t- semi-global algorithm used for suffix-prefix and prefix-suffix alignment\n" << std::endl
                              << "It prints out the statistics of two given files such as a number of sequences, maximum, minimum and average length." << std::endl
                              << "Besides that, two random sequences from the first input file are aligned and" << std::endl
                              << "the resulting alignment score and CIGAR string are printed.\n" << std::endl
                              << "The program accepts files as floating arguments." << std::endl
                              << "It supports formats: \".fasta\", \".fa\", \".fastq\", \".fq\", \".fasta.gz\", \".fa.gz\", \".fastq.gz\", \".fq.gz\"\n" << std::endl
                              << "Correct usage: blue_mapper [-m matchValue] [-s mismatchValue] [-g gapValue] [-t alignmentType] <fileame> <filename>" << std::endl
                              << "Alignment type determines which algorithm will be used. Accepted values are: \"global\", \"local\" or \"semi_global\".\n" << std::endl
                              << "It also supports options:" << std::endl
                              << "    -h --help for help menu" << std::endl
                              << "    -v --version for current version\n" << std::endl
                              << "If not provided, default values are used:" << std::endl
                              << "\t- Needleman-Wunsch algorithm" << std::endl
                              << "\t- match: 0" << std::endl
                              << "\t- mismatch: -1" << std::endl
                              << "\t- gap: -1\n" << std::endl;
                    return(0);
                case 'v':
                    std::cout << "v0.1.0" << std::endl ;
                    return(0);
                case 'm':
                    match = atoi(optarg);
                    break;
                case 's':
                    mismatch = atoi(optarg);
                    break;
                case 'g':
                    gap = atoi(optarg);
                    break;
                case 't':
                    align_type = optarg;
                    break;
                case 'k':
                    kmer_length = atoi(optarg);
                    break;
                case 'w':
                    window_length = atoi(optarg);
                    break;
                default:
                    std::cout << "The option you entered is unknown!" << std::endl;
                    exit(1);
            }
        }

        std::string type (align_type);

        if (argc != (optind + 2)) {
            std::cout << "You should've entered two files to work with. Please try again or ask for --help." << std::endl;
            exit(1);
        }

        std::vector<std::string> extensions {".fasta", ".fa", ".fastq", ".fq", ".fasta.gz", ".fa.gz", ".fastq.gz", ".fq.gz"};
        bool okay = false;

        std::string first = argv[optind];
        std::string second = argv[optind+1];


        for (size_t i=0; i<extensions.size(); i++) {
            if(first.find(extensions[i])>=0) {
                for (size_t j=0; j<extensions.size(); j++) {
                        if (second.find(extensions[j])>=0){
                            okay=true;
                            break;
                        }
                }
                if (okay) break;
            }
        }

        if (!okay){
            std::cout << "Format you entered is not compatible with our parser." << std::endl;
            exit(0);
        }

        std::vector<std::unique_ptr<InputFile>> first_object;
        if (first.find("fastq")>first.length() && first.find("fq")>first.length()){

            auto fasta_parser = bioparser::createParser<bioparser::FastaParser, InputFile>(first);
            fasta_parser->parse_objects(first_object, -1);

            std::cout << "We've parsed first file." << std::endl;
            fastaq_stat(first_object);

        } else {
            auto fastq_parser = bioparser::createParser<bioparser::FastqParser, InputFile>(first);
            uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
                while (true) {
                    auto status = fastq_parser->parse_objects(first_object, size_in_bytes);

                    if (status == false) {
                        break;
                    }
                }
            std::cout << "We've parsed first file." << std::endl;

            fastaq_stat(first_object);
        }

        std::vector<std::unique_ptr<InputFile>> second_object;
        auto fasta_parser = bioparser::createParser<bioparser::FastaParser, InputFile>(second);
        fasta_parser->parse_objects(second_object, -1);

        std::cout << "We've parsed second file." << std::endl;
        fastaq_stat(second_object);

        srand(time(NULL));
        int random_1 = rand() % first_object.size();
        int random_2 = rand() % first_object.size();

        std::string& query = first_object[random_1]->sequence;
        std::string& target = first_object[random_2]->sequence;
        std::string cigar;
        unsigned int target_begin;

        std::cout << "Alignment score: "
                  << blue::pairwise_alignment(query.c_str(), query.size(), target.c_str(), target.size(), blue::getType(type), match, mismatch, gap, cigar, target_begin) << std::endl;
        std::cout << "Cigar string: " << cigar << std::endl;
        std::cout << "Target begin: " << target_begin << std::endl;

        std::unordered_map<unsigned int, int> occurances;
        int j = 0;
        for(auto& i : first_object) {
            std::vector<std::tuple<unsigned int, unsigned int, bool>> sequenceMinimizers = blue::minimizers(i->sequence.c_str(), (i->sequence).length(), kmer_length, window_length);
            for (auto& minimizer : sequenceMinimizers)
                ++occurances[std::get<0>(minimizer)];
        }

        std::ofstream myfile;
        myfile.open ("MinimizerOccurences.csv");
        myfile << "Minimizer,Occurences\n";

        using iterator = std::unordered_map< unsigned int, int >::iterator;
        for (iterator iter = occurances.begin(); iter != occurances.end(); ++iter) {
            myfile << iter->first << "," << iter->second << '\n';
        }
        myfile.close();
        std::cout << ".CSV file with minimizer occurences is created!" << std::endl;
        return 0;
}

