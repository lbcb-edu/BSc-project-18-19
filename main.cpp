#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <memory>
#include <stdint.h>


using namespace std;

        class Example1{

            public:
                Example1(
                    const char* name, uint32_t name_length,
                    const char* sequence, uint32_t sequence_length){}

        };

        class Example2 {
            public:

            Example2(
                const char* name, uint32_t name_length,
                const char* sequence, uint32_t sequence_length,
                const char* quality, uint32_t quality_length){}
        };

int main( int argc, char* argv[] )
{

    if (argc == 2){
        int c ;
        while( ( c = getopt (argc, argv, "hv") ) != -1 )
        {
            switch(c)
            {
                case 'h':
                    cout << "You've asked for help. This is help." << endl;
                    break;
                case 'v':
                    cout << "v0.1.0" << endl ;
                    break;
                default:
                    cout << "We're starting parsing." << endl;
            }
        }

    } else if (argc == 3) {
        string extensions [8] = {".fasta", ".fa", ".fastq", ".fq", ".fasta.gz", ".fa.gz", ".fastq.gz", ".fq.gz"};
        bool okay = false;

        string first = argv[1];
        string second = argv[2];


        for (int i=0; i<8; i++) {
            if(first.find(extensions[i])>=0) {
                for (int j=0; j<8; j++) {
                        if (second.find(extensions[j])>=0){
                            okay=true;
                            break;
                        }
                }

            }
        }

        if (okay==false){
            cout << "Format you entered is not compatible with our parser." << endl;
            return 1;
        }

        if (first.find("fastq")<0 && second.find("fq")<0){

            vector<unique_ptr<Example1>> fasta_objects;
            auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Example1>(first);
            fasta_parser->parse_objects(fasta_objects, -1);
        }
        else {

            vector<unique_ptr<Example2>> fastq_objects;
            auto fastq_parser = bioparser::createParser<bioparser::FastqParser, Example2>(first);

        }

        vector<unique_ptr<Example2>> fastq_objects;
        auto fastq_parser = bioparser::createParser<bioparser::FastqParser, Example2>(second);

        uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
        while (true) {
            auto status = fastq_parser->parse_objects(fastq_objects, size_in_bytes);

            if (status == false) {
                break;
            }
        }

    }

    return 0;
}











