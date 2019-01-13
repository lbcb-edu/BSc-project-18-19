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
#include <functional>
#include <algorithm>
#include <map>
#include <set>
#include <unordered_map>
#include <cmath>
#include "bioparser/bioparser.hpp"
#include "blue_alignment.hpp"
#include "blue_minimizers.hpp"

typedef std::vector<std::tuple<unsigned int, unsigned int, bool>> uubtuple;
uubtuple findInGenome(uubtuple& sequenceMinimizers, std::unordered_map<unsigned int, uubtuple>& mapGenome);
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

std::vector<std::tuple<int, int, int, int, bool>> longestIncreasingSubSequence(uubtuple input){
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
                    len = len + 1;
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
                    len2 = len2 + 1;
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
        std::cout << "Longest increasing subsequences " << std::endl;
        int index = T[len];
        int indexPrvi = T[len];
        int indexZadnji = T[len];

        /*std::cout << "gledaj mene" << std::endl;
        while(index != -1) {
           std::cout << std::get<0>(input[index]) <<"," << std::get<1>(input[index]) << "," << std::get<2>(input[index]) << std::endl;
            index = R[index];
        }*/

        index = T[len];

        while(index != -1) {
            if (R[index] != -1 && std::abs((std::get<0>(input[index]) - std::get<0>(input[R[index]])) - (std::get<1>(input[index]) - std::get<1>(input[R[index]])) > 500))
                {
                    if(std::get<0>(input[index]) != std::get<0>(input[indexPrvi])) {
                        regions.emplace_back(std::get<0>(input[index]), std::get<0>(input[indexPrvi]), std::get<1>(input[index]), std::get<1>(input[indexPrvi]), 0);
                    }
                    indexPrvi = R[index];
                }
            indexZadnji = index;
            index = R[index];
        }

        if(index != indexPrvi) {
            regions.emplace_back(std::get<0>(input[indexZadnji]), std::get<0>(input[indexPrvi]), std::get<1>(input[indexZadnji]), std::get<1>(input[indexPrvi]), 0);
        }

        /*std::cout << "tu smo" << std::endl;
        for(auto& j : regions) {
            std::cout << "(" << std::get<0>(j) << "," << std::get<1>(j) << "," << std::get<2>(j) << "," << std::get<3>(j) << "," << std::get<4>(j) << ") " <<std::endl;
        }*/

        index = T2[len2];
        indexPrvi = T2[len2];
        indexZadnji = T2[len2];

         while(index != -1) {
            if (R2[index] != -1 && std::abs((std::get<0>(input[index]) - std::get<0>(input[R2[index]])) - (std::get<1>(input[index]) - std::get<1>(input[R2[index]])) > 500))
                {
                    if(std::get<0>(input[index]) != std::get<0>(input[indexPrvi])) {
                        regions.emplace_back(std::get<0>(input[index]), std::get<0>(input[indexPrvi]), std::get<1>(input[index]), std::get<1>(input[indexPrvi]), 1);
                    }
                    indexPrvi = R2[index];
                }
            indexZadnji = index;
            index = R2[index];
        }

        return regions;
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
    {"f", required_argument, NULL, 'f'},
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

        unsigned int kmer_length = 15;
        unsigned int window_length = 5;
        double f = 0.001;

        while ((c = getopt_long (argc, argv, "hvm:s:g:t:k:w:f:", long_options, NULL)) != -1) {
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

/*      std::unordered_map<unsigned int, int> occurences;
        int j = 0;
        for(auto& i : first_object) {
            std::vector<std::tuple<unsigned int, unsigned int, bool>> sequenceMinimizers = blue::minimizers(i->sequence.c_str(), (i->sequence).length(), kmer_length, window_length);
            for (auto& minimizer : sequenceMinimizers)
                ++occurences[std::get<0>(minimizer)];
                std::cout << ++j << std::endl;
        }

        int distinctCounter = 0;
        std::vector<std::pair<unsigned int, int>> nonDistinctMinimizers;
       	std::unordered_map<unsigned int, int>::iterator it = occurences.begin();
        while (it != occurences.end()) {
            if(it->second == 1) ++distinctCounter;
            else {
                nonDistinctMinimizers.emplace_back(it->first, it->second);
            }
            it++;
        }
        std::cout << "Number of distinct minimizers: " << distinctCounter << std::endl;

        typedef std::function<bool(std::pair<unsigned int, int>, std::pair<unsigned int, int>)> Comparator;
        Comparator comparator =
			[](std::pair<unsigned int, int> elem1, std::pair<unsigned int, int> elem2)
			{
				return elem1.second > elem2.second;
			};

        sort(nonDistinctMinimizers.begin(), nonDistinctMinimizers.end(), comparator);

        std::cout << nonDistinctMinimizers.size() << std::endl;
        int minimizer = f * nonDistinctMinimizers.size();
        std::cout << "Number of occurences of the most frequent minimizer (without top f frequent minimizers): " << nonDistinctMinimizers[minimizer].second << std::endl;*/

        uubtuple genomeMinimizers = blue::minimizers(second_object[0]->sequence.c_str(), (second_object[0]->sequence).length(), kmer_length, window_length);
        std::unordered_map<unsigned int, uubtuple> mapByValue = makeMap(genomeMinimizers);

        for(auto& i : first_object) {
            uubtuple sequenceMinimizers = blue::minimizers(i->sequence.c_str(), (i->sequence).length(), kmer_length, window_length);
            uubtuple result = findInGenome(sequenceMinimizers, mapByValue);

            sort(result.begin(), result.end(), comparator);
            std::vector<std::tuple<int, int, int, int, bool>> positions = longestIncreasingSubSequence(result);

            for(auto& position : positions){
                unsigned int querySize = std::get<1>(position)-std::get<0>(position)+1;
                unsigned int targetSize = std::get<3>(position)-std::get<2>(position)+1;

                std::string& aq = i->sequence;
                std::string& tq = (second_object[0]->sequence);

                std::string queryString = aq.substr(std::get<0>(position), querySize);
                std::string targetString = tq.substr(std::get<2>(position), targetSize);

                std::string cigar1;
                unsigned int target_begin1;

                std::cout << blue::pairwise_alignment(queryString.c_str(), querySize, targetString.c_str(), targetSize, blue::getType(type), match, mismatch, gap, cigar1, target_begin1) << std::endl;
                std::cout << cigar1 << std::endl;
                std::cout << target_begin1 << std::endl;
            }
        }
        return 0;
}

typedef std::vector<std::tuple<unsigned int, unsigned int, bool>> uubtuple;

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
