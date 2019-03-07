#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <memory>
#include <stdint.h>
#include <getopt.h>
#include <functional>
#include <algorithm>
#include <map>
#include <set>
#include <unordered_map>
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "blue_alignment.hpp"
#include "blue_minimizers.hpp"

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

typedef std::vector<std::tuple<unsigned int, unsigned int, bool>> uubtuple;
uubtuple findInGenome(uubtuple& sequenceMinimizers, std::unordered_map<unsigned int, uubtuple>& mapGenome);
void finalCountdown(std::vector<std::unique_ptr<InputFile>>& first_object, std::vector<std::unique_ptr<InputFile>>& second_object, std::unordered_map<unsigned int, uubtuple> mapByValue,
            std::string type, int thread_begin, int thread_end);
std::unordered_map<unsigned int, uubtuple> makeMap(uubtuple genome);
namespace std {
template <> struct hash<std::tuple<unsigned int, unsigned int, bool >> {
    inline size_t operator()(const std::tuple<unsigned int, unsigned int, bool > &v) const {
        std::hash<int> int_hasher;
        return int_hasher(std::get<0>(v)) ^ int_hasher(std::get<1>(v)) ^ int_hasher(std::get<2>(v));
    }
};
}


typedef std::function<bool(std::tuple<unsigned int, unsigned int, bool>, std::tuple<unsigned int, unsigned int, bool>)> Comparator;
            Comparator comparator =
                [](std::tuple<unsigned int, unsigned int, bool> elem1, std::tuple<unsigned int, unsigned int, bool> elem2)
                {
                    if(std::get<2>(elem1) == std::get<2>(elem2)) {
                        return std::get<0>(elem1) < std::get<0>(elem2);

                    } else {
                        return std::get<2>(elem1) < std::get<2>(elem2);
                    }
};

int ceilIndex(uubtuple input, int T[], int end, int s) {
        int start = 0;
        int middle;
        int len = end;
        while(start <= end) {
            middle = (start + end)/2;
            if(middle < len && std::get<1>(input[T[middle]]) < s && s <= std::get<1>(input[T[middle+1]])) {
                return middle+1;
            }else if(std::get<1>(input[T[middle]]) < s){
                start = middle+1;
            }else{
                end = middle-1;
            }
        }
        return -1;
}

std::tuple<int, int, int, int, bool> longestIncreasingSubSequence(uubtuple input){
        int T[input.size()];
        int R[input.size()];

        int T2[input.size()];
        int R2[input.size()];

        for(int i=0; i < input.size() ; i++) {
            R[i] = -1;
            R2[i] = -1;
        }

        T[0] = 0;
        int len = 0;

        T2[0] = 0;
        int len2 = 0;

        for(int i=1; i < input.size(); i++){

            if(std::get<2>(input[i]) == 0) {
                if(std::get<1>(input[T[0]]) > std::get<1>(input[i]) || std::get<0>(input[T[0]]) > std::get<0>(input[i])) { //if input[i] is less than 0th value of T then replace it there.*/
                    T[0] = i;
                }else if(std::get<1>(input[T[len]]) < std::get<1>(input[i]) &&  std::get<0>(input[T[len]]) < std::get<0>(input[i])) { //if input[i] is greater than last value of T then append it in T*/
                    len += 1;
                    T[len] = i;
                    R[T[len]] = T[len-1];

                }else{ //do a binary search to find ceiling of input[i] and put it there.
                    int index = ceilIndex(input, T, len, std::get<1>(input[i]));
                    T[index] = i;
                    R[T[index]] = T[index-1];
                }

            } else {
                 if(std::get<1>(input[T2[0]]) > std::get<1>(input[i]) || std::get<0>(input[T2[0]]) > std::get<0>(input[i])) { //if input[i] is less than 0th value of T then replace it there.
                    T2[0] = i;
                }else if(std::get<1>(input[T2[len2]]) < std::get<1>(input[i]) &&  std::get<0>(input[T2[len2]]) < std::get<0>(input[i])) { //if input[i] is greater than last value of T then append it in T
                    len2 += 1;
                    T2[len2] = i;
                    R2[T2[len2]] = T2[len2-1];

                }else{ //do a binary search to find ceiling of input[i] and put it there.
                    int index = ceilIndex(input, T2, len2, std::get<1>(input[i]));
                    T2[index] = i;
                    R2[T2[index]] = T2[index-1];
                }
            }
        }

        std::vector<std::tuple<int, int, int, int, bool>> regions;
        int index = T[len];
        int pocetniIndex = T[len];
        int krajnjiIndex = T[len];

        while(index != -1){
            int absV = (std::get<1>(input[index]) - std::get<0>(input[index])) - (std::get<1>(input[R[index]]) - std::get<0>(input[R[index]]));
            absV = absV < 0 ? (-1)*absV : absV;
            if (absV < 500) {
                krajnjiIndex = R[index];
            } else {
                if (std::get<0>(input[pocetniIndex]) - std::get<0>(input[krajnjiIndex]) > 4 && std::get<1>(input[pocetniIndex]) - std::get<1>(input[krajnjiIndex]) > 4)
                    regions.emplace_back(std::get<0>(input[krajnjiIndex]), std::get<0>(input[pocetniIndex]), std::get<1>(input[krajnjiIndex]), std::get<1>(input[pocetniIndex]), 0);
                pocetniIndex = R[index];
                krajnjiIndex = R[index];
            }
            index = R[index];
        }

        if (std::get<0>(input[pocetniIndex]) - std::get<0>(input[krajnjiIndex]) > 4 && std::get<1>(input[pocetniIndex]) - std::get<1>(input[krajnjiIndex]) > 4)
            regions.emplace_back(std::get<0>(input[krajnjiIndex]), std::get<0>(input[pocetniIndex]), std::get<1>(input[krajnjiIndex]), std::get<1>(input[pocetniIndex]), 0);


        index = T2[len];
        pocetniIndex = T2[len];
        krajnjiIndex = T2[len];

        while(index != -1){
            int absV = (std::get<1>(input[index]) - std::get<0>(input[index])) - (std::get<1>(input[R2[index]]) - std::get<0>(input[R2[index]]));
            absV = absV < 0 ? (-1)*absV : absV;
            if (absV < 500) {
                krajnjiIndex = R2[index];
            } else {
                if (std::get<0>(input[pocetniIndex]) - std::get<0>(input[krajnjiIndex]) > 4 && std::get<1>(input[pocetniIndex]) - std::get<1>(input[krajnjiIndex]) > 4)
                    regions.emplace_back(std::get<0>(input[krajnjiIndex]), std::get<0>(input[pocetniIndex]), std::get<1>(input[krajnjiIndex]), std::get<1>(input[pocetniIndex]), 1);
                pocetniIndex = R2[index];
                krajnjiIndex = R2[index];
            }
            index = R2[index];
        }

        if (std::get<0>(input[pocetniIndex]) - std::get<0>(input[krajnjiIndex]) > 4 && std::get<1>(input[pocetniIndex]) - std::get<1>(input[krajnjiIndex]) > 4)
            regions.emplace_back(std::get<0>(input[krajnjiIndex]), std::get<0>(input[pocetniIndex]), std::get<1>(input[krajnjiIndex]), std::get<1>(input[pocetniIndex]), 1);

        int maxLength = 0;
        std::tuple<int, int, int, int, bool> maxTuple;
        for(auto &i : regions) {
            if(std::get<1>(i) - std::get<0>(i) > 4 && std::get<3>(i) - std::get<2>(i) > 4) {
                if(std::get<3>(i) - std::get<2>(i) > maxLength){
                    maxLength = std::get<3>(i) - std::get<2>(i);
                    maxTuple = i;
                }
            } else
                continue;
        }

        return maxTuple;
}


static struct option long_options[] = {
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    {"match", required_argument, NULL, 'm' },
    {"mismatch", required_argument, NULL, 's'},
    {"gap", required_argument, NULL, 'g'},
    {"type", required_argument, NULL, 't'},
    {"kmer_length", required_argument, NULL, 'k'},
    {"window_length", required_argument, NULL, 'w'},
    {"f", required_argument, NULL, 'f'},
    {"cigar", required_argument, NULL, 'c'},
    {"paralelization", required_argument, NULL, 'p'},
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

int c;
int match = 0;
int mismatch = -1;
int gap = -1;
char* align_type = "global";

unsigned int kmer_length = 15;
unsigned int window_length = 5;
double f = 0.001;

bool cCigar = false;
int paralelization = 1;

int main (int argc, char* argv[]) {


        while ((c = getopt_long (argc, argv, "hvm:s:g:t:k:w:f:cp:", long_options, NULL)) != -1) {
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
                case 'f':
                    f = atof(optarg);
                    break;
                case 'c':
                    cCigar = true;
                    break;
                case 'p':
                    paralelization = atoi(optarg);
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

        uubtuple genomeMinimizers = blue::minimizers(second_object[0]->sequence.c_str(), (second_object[0]->sequence).length(), kmer_length, window_length);
        std::unordered_map<unsigned int, uubtuple> mapByValue = makeMap(genomeMinimizers);

        std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(paralelization);
        std::vector<std::future<void>> thread_futures;

        int each_thread_work = first_object.size() / paralelization;

        int thread_begin = 0;
        int thread_end = each_thread_work + (first_object.size() % paralelization);

        for (int i = 0; i < paralelization; ++i) {
            thread_futures.emplace_back(thread_pool->submit_task(finalCountdown, std::ref(first_object), std::ref(second_object), mapByValue, type, thread_begin, thread_end));
            thread_begin = thread_end-1;
            thread_end = thread_begin + each_thread_work;
        }

        for (auto& it: thread_futures) {
           it.wait();
        }

        /*uubtuple proba{ std::make_tuple(1, 2, false), std::make_tuple(2, 15, false), std::make_tuple(3000, 5000, false), std::make_tuple(3000, 5400, false), std::make_tuple(3500, 5900, false) };

        std::tuple<int, int, int, int, bool> result = longestIncreasingSubSequence(proba);

        std::cout << std::get<0>(result) << " " << std::get<1>(result) << " " << std::get<2>(result) << " " << std::get<3>(result) << std::endl;*/


        return 0;
}

std::unordered_map<unsigned int, uubtuple> makeMap(uubtuple genome) {
    std::unordered_map<unsigned int, uubtuple> mapa;
    for(auto &i : genome) {
        std::unordered_map<unsigned int, uubtuple>::const_iterator got = mapa.find(std::get<0>(i));
        if ( got == mapa.end() ) {
            uubtuple novi;
            novi.push_back(i);
            mapa[std::get<0>(i)] = novi;

        } else {
            uubtuple novi2 = got->second;
            novi2.push_back(i);
            mapa[std::get<0>(i)] = novi2;
        }
    }
    return mapa;
}

uubtuple findInGenome(uubtuple& sequenceMinimizers, std::unordered_map<unsigned int, uubtuple>& mapGenome) {
    uubtuple result;
    for(auto& i : sequenceMinimizers) {
        std::unordered_map<unsigned int, uubtuple>::const_iterator got = mapGenome.find(std::get<0>(i));
        if(got == mapGenome.end()) {
            continue;

        } else {
            for(auto& j : got->second) {
                result.emplace_back(std::get<1>(i), std::get<1>(j), std::get<2>(i) ^ std::get<2>(j));
            }
        }
    }
    return result;
}

void finalCountdown(std::vector<std::unique_ptr<InputFile>>& first_object, std::vector<std::unique_ptr<InputFile>>& second_object, std::unordered_map<unsigned int, uubtuple> mapByValue,
                        std::string type, int thread_begin, int thread_end) {
    int j=0;
    for(int z = thread_begin; z < thread_end; z++) {
            ++j;
            auto &i = first_object[z];
            std::cout << j << std::endl;
            uubtuple sequenceMinimizers = blue::minimizers(i->sequence.c_str(), (i->sequence).length(), kmer_length, window_length);
            //continue;

            uubtuple result = findInGenome(sequenceMinimizers, mapByValue);

            if (result.size() == 0) {
                continue;
            }

            sort(result.begin(), result.end(), comparator);

            std::cout << "tu" << std::endl;

            std::tuple<int, int, int, int, bool> position = longestIncreasingSubSequence(result);

            std::cout << std::get<0>(position) << ", " << std::get<1>(position) << ", " << "blabla" << std::endl;

            if(std::get<0>(position) == 0 && std::get<1>(position) == 0) continue;

            unsigned int querySize = std::get<1>(position)-std::get<0>(position)+1;
            unsigned int targetSize = std::get<3>(position)-std::get<2>(position)+1;

            std::string& aq = i->sequence;
            std::string& tq = (second_object[0]->sequence);

            std::string queryString = aq.substr(std::get<0>(position), querySize);
            std::string targetString = tq.substr(std::get<2>(position), targetSize);

            std::string cigar1;
            unsigned int target_begin1;

            blue::pairwise_alignment(queryString.c_str(), querySize, targetString.c_str(), targetSize, blue::getType(type), match, mismatch, gap, cigar1, target_begin1);

            //PAF FORMAT
            std::string paf = i->name + '\t' + std::to_string((i->sequence).size()) + '\t' + std::to_string(std::get<0>(position)) + '\t' + std::to_string(std::get<1>(position)) + '\t';
            paf += std::get<4>(position) == 0 ? '+' : '-';
            paf = paf + '\t' + second_object[0]->name + '\t' + std::to_string((second_object[0]->sequence).size()) + '\t' + std::to_string(std::get<2>(position));
            paf = paf + '\t' + std::to_string(std::get<3>(position)) +'\t';

            std::string number;
            int noOfMatches = 0;
            int blockLength = 0;

            for(char c : cigar1) {
                if(isdigit(c)) {
                    number += c;
                } else if(c == '=') {
                    noOfMatches += stoi(number);
                    blockLength += stoi(number);
                    number = "";
                } else {
                    blockLength += stoi(number);
                    number = "";
                }
            }

            paf += std::to_string(noOfMatches) + '\t' + std::to_string(blockLength) + '\t';
            paf += (i->quality).empty() ? "255\t" : i->quality;

            if(cCigar) {
                paf += "cg:Z:" + cigar1;
            }

            std::cout << paf << std::endl;
        }
}
