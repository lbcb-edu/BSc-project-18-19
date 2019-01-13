#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <getopt.h>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <unordered_map>
#include <fstream>
#include <cmath>
#include <tr1/memory>
#include <iterator>

#include <bioparser/bioparser.hpp>
#include "pink_alignment.hpp"
#include "pink_minimizers.hpp"

static struct option options[] = {
    {"help",          no_argument,       0, 'h'},
    {"version",       no_argument,       0, 'v'},
    {"global",        no_argument,       0, 'G'},
    {"semi_global",   no_argument,       0, 'S'},
    {"local",         no_argument,       0, 'L'},
    {"match",         required_argument, 0, 'm'},
    {"mismatch",      required_argument, 0, 's'},
    {"gap",           required_argument, 0, 'g'},
    {"k",             required_argument, 0, 'k'},
    {"window_length", required_argument, 0, 'w'},
    {"f",             required_argument, 0, 'f'},
    {"cigar",         no_argument,       0, 'c'},
    {"thread",        required_argument, 0, 't'},
    {NULL,            no_argument,       0, 0}
};

class Fast {
public:
    std::string name;
    std::string sequence;
    std::string quality;

    Fast(
        const char *name, uint32_t name_length,
        const char *sequence, uint32_t sequence_length) :
        name{std::string(name, name_length)},
        sequence{std::string(sequence, sequence_length)} {}

    Fast(
        const char *name, uint32_t name_length,
        const char *sequence, uint32_t sequence_length,
        const char *quality, uint32_t quality_length) :
        name{std::string(name, name_length)},
        sequence{std::string(sequence, sequence_length)},
        quality{std::string(quality, quality_length)} {}
};

std::vector<std::unique_ptr<Fast>> parse_fasta(std::string fastaFile) {
    std::vector<std::unique_ptr<Fast>> fasta_objects;

    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Fast>(fastaFile);
    fasta_parser->parse_objects(fasta_objects, -1);

    return fasta_objects;
}

std::vector<std::unique_ptr<Fast>> parse_fastq(std::string fastqFile) {
    std::vector<std::unique_ptr<Fast>> fastq_objects;

    auto fastq_parser = bioparser::createParser<bioparser::FastqParser, Fast>(fastqFile);
    uint64_t size_in_bytes = 500 * 1024 * 1024;  // 500 MB

    while (true) {
        auto status = fastq_parser->parse_objects(fastq_objects, size_in_bytes);
        if (status == false) {
            break;
        }
    }

    return fastq_objects;
}

void print_stats(const std::vector<std::unique_ptr<Fast>> &fast_objects) {
    unsigned numOfSeq = fast_objects.size();
    unsigned sum = 0;
    float average;
    unsigned min = (fast_objects[0]->sequence).length();
    unsigned max = min;

    for (unsigned i = 0; i < numOfSeq; i++) {
        sum += (fast_objects[i]->sequence).length();

        if ((fast_objects[i]->sequence).length() < min) {
            min = (fast_objects[i]->sequence).length();
        }
        if ((fast_objects[i]->sequence).length() > max) {
            max = (fast_objects[i]->sequence).length();
        }
    }
    average = (float) sum / numOfSeq;

    std::cerr << "Number of sequences: " << numOfSeq << std::endl;
    std::cerr << "Average length: " << average << std::endl;
    std::cerr << "Minimal length: " << min << std::endl;
    std::cerr << "Maximal length: " << max << std::endl;
}

bool check_extension(std::string arg, char ext_flag) {
    std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);

    std::vector<std::string> ext;
    if (ext_flag == 'a') {
        ext = {".fa", ".fasta", ".fa.gz", ".fasta.gz"};
    } else {
        ext = {".fq", ".fastq", ".fq.gz", ".fastq.gz"};
    }

    for (unsigned i = 0; i < ext.size(); i++) {
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

void printError() {
    std::cerr << "Wrong input. Use \"-h\" or \"--help\" for help." << std::endl;
    exit(1);
}

int checkInput(char *optarg) {
    if (atoi(optarg) == 0 && strcmp(optarg, "0") != 0) {
        printError();
    }
    return atoi(optarg);
}

void help() {
    printf(
            "usage: pink_mapper [options ...] <fragments> <genome> \n"
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
            "           prints the usage\n"
            "       -G, --global\n"
            "           set to global alignment, default alignment\n"
            "       -S, --semi_global\n"
            "           set to semi-global alignment\n"
            "       -L, --local\n"
            "           set to local alignment\n"
            "       -m, --match\n"
            "           input match cost, default: 2\n"
            "       -s, --mismatch\n"
            "           input mismatch cost, default: -1\n"
            "       -g, --gap\n"
            "           input insertion/deletion cost, default: -2\n"
            "       -k\n"
            "           input k, default: 15\n"
            "       -w\n"
            "           input window length, default: 5\n"
            "       -f\n"
            "           input percentage of minimizers which are not taken in account\n"
            "		-c, --cigar\n"
            "			include CIGAR strings of the alignment\n"
            "		-t, --thread\n"
            "			input number of threads\n"
    );
}

void version() {
    printf("v0.1.0\n");
}

template<typename E>
struct Node {
    E value;
    std::tr1::shared_ptr<Node<E>> pointer;
};

template<class E>
struct node_ptr_less {
    bool operator()(const std::tr1::shared_ptr<Node<E> > &node1,
                    const std::tr1::shared_ptr<Node<E> > &node2) const {
        return node1->value < node2->value;
    }
};

template<typename E>
std::vector<E> lis(const std::vector<E> &n) {
    typedef std::tr1::shared_ptr<Node<E> > NodePtr;

    std::vector<NodePtr> pileTops;
    // sort into piles
    for (typename std::vector<E>::const_iterator it = n.begin(); it != n.end(); it++) {
        NodePtr node(new Node<E>());
        node->value = *it;
        typename std::vector<NodePtr>::iterator j =
                std::lower_bound(pileTops.begin(), pileTops.end(), node, node_ptr_less<E>());
        if (j != pileTops.begin())
            node->pointer = *(j - 1);
        if (j != pileTops.end())
            *j = node;
        else
            pileTops.push_back(node);
    }
    // extract LIS from piles
    std::vector<E> result;
    for (NodePtr node = pileTops.back(); node != NULL; node = node->pointer)
        result.push_back(node->value);
    std::reverse(result.begin(), result.end());
    return result;
}

bool sortbysec(const std::pair<int, int> &a,
               const std::pair<int, int> &b) {
    return (a.second < b.second);
}


int main(int argc, char *argv[]) {
    char optchr;
    srand(time(NULL));

    pink::AlignmentType type = pink::global;
    int match = 2;
    int mismatch = -1;
    int gap = -2;
    std::string cigar;
    unsigned int target_begin = 0;

    unsigned int k = 15;
    unsigned int window_length = 5;
    float f = 0.001;

    //bool c = false;
    //int thread = 1;

    while ((optchr = getopt_long(argc, argv, "hvGSLm:s:g:k:w:f:ct:", options, NULL)) != -1) {
        switch (optchr) {
            case 'h':
                help();
                return 0;
            case 'v':
                version();
                return 0;
            case 'G':
                type = pink::global;
                break;
            case 'S':
                type = pink::semi_global;
                break;
            case 'L':
                type = pink::local;
                break;
            case 'm':
                match = checkInput(optarg);
                break;
            case 's':
                mismatch = checkInput(optarg);
                break;
            case 'g':
                gap = checkInput(optarg);
                break;
            case 'k':
                k = checkInput(optarg);
                break;
            case 'w':
                window_length = checkInput(optarg);
                break;
            case 'f':
                f = atof(optarg);
                break;
            case 'c':
                //c = true;
                break;
            case 't':
                //thread = checkInput(optarg);
                break;
            default:
                printError();
        }
    }

    if (argc - optind != 2) {
        printError();
    }

    if ((check_extension(argv[optind], 'a') || check_extension(argv[optind], 'q')) &&
        check_extension(argv[optind + 1], 'a')) {
        std::vector<std::unique_ptr<Fast>> fast_objects1;
        std::vector<std::unique_ptr<Fast>> fast_objects2;

        if (check_extension(argv[optind], 'a')) {
            fast_objects1 = parse_fasta(argv[optind]);
        } else {
            fast_objects1 = parse_fastq(argv[optind]);
        }

        fast_objects2 = parse_fasta(argv[optind + 1]);

        std::cerr << "~FIRST FILE~" << std::endl;
        print_stats(fast_objects1);
        std::cerr << "\n" << "~SECOND FILE~" << std::endl;
        print_stats(fast_objects2);


        int query = rand() % fast_objects1.size();
        int target = rand() % fast_objects1.size();

        const char *q = (fast_objects1[query]->sequence).c_str();
        const char *t = (fast_objects1[target]->sequence).c_str();
        unsigned int q_len = (fast_objects1[query]->sequence).length();
        unsigned int t_len = (fast_objects1[target]->sequence).length();

        int cost = pink::pairwise_alignment(q, q_len, t, t_len, type, match, mismatch, gap);
        std::cout << "\nFinal cost: " << cost << std::endl;

        pink::pairwise_alignment(q, q_len, t, t_len, type, match, mismatch, gap, cigar, target_begin);
        cigar = std::string(cigar.rbegin(), cigar.rend());
        std::cout << "\nCigar: " << cigar << std::endl;



    //create minimizer index from the reference genome
        std::cout << "\nPlease wait until the minimizers are processed..." << std::endl;

        std::vector<std::tuple<unsigned int, unsigned int, bool>> t_minimizer_vector;
        std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> t_minimizer_index;
        std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>>::iterator t_it;

        t_minimizer_vector = pink::minimizers((fast_objects2.front()->sequence).c_str(),
                                              (fast_objects2.front()->sequence).length(), k, window_length);

        for (auto const &minimizer: t_minimizer_vector) {
            t_it = t_minimizer_index.find(std::get<0>(minimizer));

            if (t_it == t_minimizer_index.end()) {
                std::vector<std::pair<unsigned int, bool>> positions;
                positions.emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
                t_minimizer_index.insert(
                        std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>>(std::get<0>(minimizer),
                                                                                            positions));

            } else {
                (t_it->second).emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
            }
        }

        std::vector<std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>>> temp_index(t_minimizer_index.begin(),
                                                                                               t_minimizer_index.end());

        auto comparator = [](std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>> a,
                             std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>> b) {
            return (a.second).size() > (b.second).size();
        };

        std::sort(temp_index.begin(), temp_index.end(), comparator);

        unsigned int x = std::round(f * t_minimizer_index.size());
        temp_index.erase(temp_index.begin(), temp_index.begin() + x);

        std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> t_index(temp_index.begin(),
                                                                                                        temp_index.end());
        std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>>::iterator it;

    //create minimizer index for each fragment
        std::vector<std::tuple<unsigned int, unsigned int, bool>> q_minimizer_vector;
        std::vector<std::pair<unsigned int, unsigned int>> same;
        std::vector<std::pair<unsigned int, unsigned int>> different;

        type = pink::local;
        std::string t_v2 = fast_objects2[0]->sequence;

        for (auto const &query : fast_objects1) {

            q_minimizer_vector = pink::minimizers((query->sequence).c_str(),
                                                  (query->sequence).length(), k, window_length);

            for (auto const &minimizer: q_minimizer_vector) {
                it = t_index.find(std::get<0>(minimizer));

                if (it != t_index.end()) {
                    for (auto i : it->second) {

                        if ((std::get<2>(minimizer) && std::get<1>(i)) ||
                            (!std::get<2>(minimizer) && !std::get<1>(i))) {
                            same.emplace_back(std::make_pair(std::get<0>(i), std::get<1>(minimizer)));

                        } else {
                            different.emplace_back(std::make_pair(std::get<0>(i), std::get<1>(minimizer)));
                        }
                    }
                }
            }

            sort(same.begin(), same.end(), sortbysec);
            sort(different.begin(), different.end(), sortbysec);

        //LIS
            std::vector<std::pair<unsigned int, unsigned int>> s_locations = lis(same);
            std::vector<std::pair<unsigned int, unsigned int>> d_locations = lis(different);

//            for (auto s : same) {
//                std::cout << std::get<0>(s) << ", " << std::get<1>(s) << std::endl;
//            }
//
//            for (auto s : s_locations){
//                std::cout << std::get<0>(s) << ", " << std::get<1>(s) << std::endl;
//            }

            unsigned int q_begin;
            unsigned int q_end;
            unsigned int t_begin;
            unsigned int t_end;

            if (s_locations.size() >= d_locations.size()) {
                q_begin = std::get<1>(s_locations.at(0));
                q_end = std::get<1>(s_locations.at(s_locations.size() - 1));
                t_begin = std::get<0>(s_locations.at(0));
                t_end = std::get<0>(s_locations.at(s_locations.size() - 1));
            } else {
                q_begin = std::get<1>(d_locations.at(0));
                q_end = std::get<1>(d_locations.at(s_locations.size() - 1));
                t_begin = std::get<0>(d_locations.at(0));
                t_end = std::get<0>(d_locations.at(s_locations.size() - 1));
            }

            std::string q_v2 = query->sequence;

            q_len = q_end - q_begin + 1;
            t_len = t_end - t_begin + 1;


//            std::cout << "q_begin: " << q_begin << std::endl;
//            std::cout << "q_end: " << q_end << std::endl;
//            std::cout << "q_len: " << q_len << std::endl;
//            std::cout << "q_v2_len:" << q_v2.length() << std::endl;

            std::string q_v3 = q_v2.substr(q_begin, q_len);
//            std::cout << "q_v3: " << q_v3 << std::endl;

//            std::cout << "t_begin: " << t_begin << std::endl;
//            std::cout << "t_end: " << t_end << std::endl;
//            std::cout << "t_len: " << t_len << std::endl;
//            std::cout << "t_v2_len:" << t_v2.length() << std::endl;

            std::string t_v3 = t_v2.substr(t_begin, t_len);
//            std::cout << "t_v3: " << t_v3 << std::endl;


            pink::pairwise_alignment(q_v3.c_str(), q_len, t_v3.c_str(), t_len, type, match, mismatch, gap, cigar, target_begin);
            cigar = std::string(cigar.rbegin(), cigar.rend());
            std::cout << "\nCigar: " << cigar << std::endl;


            same.clear();
            different.clear();
        }

    } else {
        printError();
    }

    return 0;
}  
