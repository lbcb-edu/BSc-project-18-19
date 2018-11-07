#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <memory>
#include <stdint.h>
#include <getopt.h>
#include "bioparser/bioparser.hpp"
#include "blue_alignment.hpp"

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
    std::cout << "Number of sequences is: "<< fq_objects.size() << std::endl;
};

int main (int argc, char* argv[]) {

        int c;
        int match = 1;
        int mismatch = 1;
        int gap = -1;
        char* align_type = "global";

        while ((c = getopt_long (argc, argv, "hvm:s:g:t:", long_options, NULL)) != -1) {
            switch(c) {
                case 'h':
                    std::cout << "You've asked for help. " << std::endl << "This program accepts two files as floating arguments." << std::endl
                    << "It prints out the statistics such as number of sequnces, maximum, minimun and average length." << std::endl
                    << "It supports formats: \".fasta\", \".fa\", \".fastq\", \".fq\", \".fasta.gz\", \".fa.gz\", \".fastq.gz\", \".fq.gz\"" << std::endl
                    << "correct usage: blue_mapper <fileame> <filename>" << std::endl
                    << "It also supports options:" << std::endl
                    << "    -h --help for help menu" << std::endl
                    << "    -v --version for current version" << std::endl;
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

        std::string query = first_object.at(random_1)->sequence;
        std::string target = first_object.at(random_2)->sequence;
        unsigned int query_length = first_object.at(random_1)->sequence.length();
        unsigned int target_length = first_object.at(random_2)->sequence.length();
        std::string cigar;
        unsigned int target_begin;

        std::cout << blue::pairwise_alignment(query.c_str(), query_length, target.c_str(), target_length, blue::getType(type), match, mismatch, gap, cigar, target_begin) << std::endl;
        std::cout << cigar;

    return 0;
}

