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

class Example1{

    public:
    const char* name;
    uint32_t name_length;
    const char* sequence;
    uint32_t sequence_length;

    Example1(
        const char* name, uint32_t name_length,
        const char* sequence, uint32_t sequence_length
    ):  name{name}, name_length{name_length}, sequence{sequence}, sequence_length{sequence_length} {
        }

};

class Example2 {
    public:
    const char* name;
    uint32_t name_length;
    const char* sequence;
    uint32_t sequence_length;
    const char* quality;
    uint32_t quality_length;

    Example2(
        const char* name, uint32_t name_length,
        const char* sequence, uint32_t sequence_length,
        const char* quality, uint32_t quality_length
    ):  name{name}, name_length{name_length}, sequence{sequence},
        sequence_length{sequence_length}, quality{quality}, quality_length{quality_length} {
        }
};

static struct option long_options[] = {
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    {NULL, no_argument, NULL, 0}
};

template<typename T>
void fastaq_stat(std::vector<std::unique_ptr<T>>& fq_objects) {

    uint32_t max=0, min=INT_MAX;
    uint64_t sum=0;

    for (auto& i : fq_objects) {
        if (i->sequence_length>max){
            max=i->sequence_length;
        }

        if (i->sequence_length<min) {
            min=i->sequence_length;
        }

        sum+=i->sequence_length;
    }

    std::cout << "Maximum length of a sequence is: " << max << std::endl;
    std::cout << "Minimum length of a sequence is: " << min << std::endl;
    std::cout << "Average length of a sequence is: " << sum/(double)fq_objects.size() << std::endl;
    std::cout << "Number of sequences is: "<< fq_objects.size() << std::endl;
};



int main (int argc, char* argv[]) {

    if (argc == 2) {
        int c;
        while ((c = getopt_long (argc, argv, "hv", long_options, NULL)) != -1) {
            switch(c) {
                case 'h':
                    std::cout << "You've asked for help. This is help." << std::endl;
                    break;
                case 'v':
                    std::cout << "v0.1.0" << std::endl ;
                    break;
                default:
                    std::cout << "The option you entered is unknown!" << std::endl;
            }
        }

    } else if (argc == 3) {
        std::vector<std::string> extensions {".fasta", ".fa", ".fastq", ".fq", ".fasta.gz", ".fa.gz", ".fastq.gz", ".fq.gz"};
        bool okay = false;

        std::string first = argv[1];
        std::string second = argv[2];


        for (size_t i=0; i<extensions.size(); i++) {

            if(first.find(extensions[i])>=0) {
                for (size_t j=0; j<extensions.size(); j++) {
                        if (second.find(extensions[j])>=0){
                            okay=true;
                            break;
                        }
                }
            if (okay) {
                break;
            }
            }
        }

        if (okay==false){
            std::cout << "Format you entered is not compatible with our parser." << std::endl;
            return 1;
        }

        if (first.find("fastq")>first.length() && first.find("fq")>first.length()){

            std::vector<std::unique_ptr<Example1>> first_object;
            auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Example1>(first);
            fasta_parser->parse_objects(first_object, -1);

            std::cout << "We've parsed first file." << std::endl;
            fastaq_stat(first_object);

        } else {
            std::vector<std::unique_ptr<Example2>> first_object;
            auto fastq_parser = bioparser::createParser<bioparser::FastqParser, Example2>(first);

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

        std::vector<std::unique_ptr<Example1>> second_object;
        auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Example1>(second);
        fasta_parser->parse_objects(second_object, -1);

        std::cout << "We've parsed second file." << std::endl;

        fastaq_stat(second_object);
    }

    return 0;
}











