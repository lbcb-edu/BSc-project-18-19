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
#include <memory>
#include <iterator>

#include <bioparser/bioparser.hpp>
#include "pink_alignment.hpp"
#include "pink_minimizers.hpp"

#define OFFSET 500

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
    std::shared_ptr<Node<E>> pointer;
};

template<class E>
struct node_ptr_less {
    bool operator()(const std::shared_ptr<Node<E> > &node1,
                    const std::shared_ptr<Node<E> > &node2) const {
        return node1->value < node2->value;
    }
};

template<typename E>
std::vector<E> lis(const std::vector<E> &n) {
    typedef std::shared_ptr<Node<E>> NodePtr;

//    std::cout << "n.size() " << n.size() << std::endl;

//    std::cout << "lis " << std::endl;

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

//    std::cout << "for lis " << std::endl;

//    std::cout << pileTops.size() << std::endl;

    // extract LIS from piles
    std::vector<E> result;
    for (NodePtr node = pileTops.back(); node != NULL; node = node->pointer)
        result.push_back(node->value);
    std::reverse(result.begin(), result.end());

//    std::cout << "end lis " << std::endl;

    return result;
}

auto comparator = [](std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>> a,
                     std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>> b) {
    return (a.second).size() > (b.second).size();
};

bool sortbysec(const std::pair<int, int> &a,
               const std::pair<int, int> &b) {
    return (a.second < b.second);
}

//create minimizer index from the reference genome
void createTargetIndex(const char *target, unsigned int t_len, unsigned int k, unsigned int w, float f,
                       std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> *t_index) {
    std::cout << "\nCreating minimizer index from the reference genome..." << std::endl;

    std::vector<std::tuple<unsigned int, unsigned int, bool>> t_minimizer_vector;
    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> t_minimizer_index;
    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>>::iterator t_it;

    t_minimizer_vector = pink::minimizers(target, t_len, k, w);

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

    std::vector<std::pair<unsigned int, std::vector<std::pair<unsigned int, bool>>>> temp_index(
            t_minimizer_index.begin(),
            t_minimizer_index.end());

    std::sort(temp_index.begin(), temp_index.end(), comparator);

    unsigned int x = std::round(f * t_minimizer_index.size());
    temp_index.erase(temp_index.begin(), temp_index.begin() + x);

    for (auto temp: temp_index) {
        (*t_index).insert(temp);
    }
}

void
matchSequences(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> *regions,
               std::vector<std::pair<unsigned int, unsigned int>> *same, unsigned int s_len,
               std::vector<std::pair<unsigned int, unsigned int>> *different, unsigned int d_len) {

    sort((*same).begin(), (*same).end(), sortbysec);
    sort((*different).begin(), (*different).end(), sortbysec);

    //LIS
    std::vector<std::pair<unsigned int, unsigned int>> s_locations;
    std::vector<std::pair<unsigned int, unsigned int>> d_locations;

    unsigned int s_locations_len = 0;
    unsigned int d_locations_len = 0;


    if (s_len != 0) {
        s_locations = lis(*same);
        s_locations_len = s_locations.size();
    }

    if (d_len != 0) {
        d_locations = lis(*different);
        d_locations_len = d_locations.size();
    }

    bool s = true;

    if (s_locations_len < d_locations_len) {
        s_locations = d_locations;
        s = false;
    }

    std::cout << "s_locations!!!" << std::endl;
    for (auto sl : s_locations) {
        std::cout << std::get<0>(sl) << ", " << std::get<1>(sl) << std::endl;
    }

    std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool> region = std::make_tuple(std::get<1>(s_locations[0]), 0, std::get<0>(s_locations[0]), 0, s);

    for (unsigned i = 1; i < s_locations.size(); i++) {
        int tmp = std::get<0>(s_locations[i]) - std::get<0>(s_locations[i-1]);

        if (tmp > OFFSET) {
            if (std::get<3>(region) != 0) {
                (*regions).emplace_back(region);
            }
            region = std::make_tuple(std::get<1>(s_locations[i]), 0, std::get<0>(s_locations[i]), 0, s);

        } else {
            std::get<3>(region) += tmp;
            std::get<1>(region) += std::get<1>(s_locations[i]) - std::get<1>(s_locations[i-1]);
        }
    }
}

void
printPAF(const char *q_name, unsigned int q_len, const char *t_name, unsigned int t_len, unsigned int k,
         std::string cigar, bool c, bool s) {

    std::string pafFormat = std::string(q_name) + '\n' + std::to_string(q_len) + '\n' + '0' + '\n' +
                            std::to_string(q_len - k) + '\n';
    if (s) {
        pafFormat += "+\n";
    } else {
        pafFormat += "-\n";
    }

    pafFormat += std::string(t_name) + '\n' + std::to_string(t_len) + '\n' + '0' + '\n' +
                 std::to_string(t_len - k) + '\n';

    int numberOfMatches = 0;
    int blockLength = cigar.size();

    for (char c: cigar) {
        if (c == '=')
            numberOfMatches++;
    }

    pafFormat += std::to_string(numberOfMatches) + '\n' + std::to_string(blockLength) + '\n' + "255" + '\n';

    if (c)
        pafFormat += cigar;

    std::cout << pafFormat << std::endl;
}

//create minimizer index for each fragment
void createQueryIndex(const std::vector<std::unique_ptr<Fast>> &fast_objects1,
                      const std::vector<std::unique_ptr<Fast>> &fast_objects2,
                      unsigned int k, unsigned int w, float f, int match, int mismatch, int gap, bool c) {
    int i = 0;
    auto target = fast_objects2.front()->sequence;

    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> t_index;
    createTargetIndex(target.c_str(), target.length(), k, w, f, &t_index);
    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>>::iterator it;

    for (auto const &query_object : fast_objects1) {
        std::cout << "/nNo. of query sequence : " << i << std::endl;
        auto query = query_object->sequence;

        std::vector<std::tuple<unsigned int, unsigned int, bool>> q_minimizer_vector = pink::minimizers(query.c_str(),
                                                                                                        query.length(),
                                                                                                        k, w);
        std::vector<std::pair<unsigned int, unsigned int>> same;
        std::vector<std::pair<unsigned int, unsigned int>> different;

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

        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, bool>> regions;
        matchSequences(&regions, &same, same.size(), &different, different.size());

        same.clear();
        same.shrink_to_fit();
        different.clear();
        different.shrink_to_fit();
        q_minimizer_vector.clear();
        q_minimizer_vector.shrink_to_fit();

        // std::string cigar = "";
//        pink::AlignmentType type = pink::local;
        std::string cigar2;
//        unsigned int target_begin = 0;

        for (auto const &region : regions) {
            std::cout << "region!!!" << std::endl;
            std::cout << "q_len: " << std::get<1>(region) << std::endl;
            std::cout << "t_len: " << std::get<3>(region) << std::endl;

            std::string q_sub = query.substr(std::get<0>(region), std::get<1>(region));
            std::string t_sub = target.substr(std::get<2>(region), std::get<3>(region));

            std:: cout << "std::get<0>(region) " << std::get<0>(region) << std::endl;
            std:: cout << "std::get<1>(region) " << std::get<1>(region) << std::endl;
            std:: cout << "std::get<2>(region) " << std::get<2>(region) << std::endl;
            std:: cout << "std::get<3>(region) " << std::get<3>(region) << std::endl;

            std::string cigar;

            std::cout << "q_sub.c_str() " << q_sub.c_str() << std::endl;
            std::cout << "t_sub.c_str() " << t_sub.c_str() << std::endl;

            //pink::pairwise_alignment(q_sub.c_str(), q_sub.size(), t_sub.c_str(), t_sub.size(), type, match, mismatch, gap, cigar, target_begin);
            //cigar = std::string(cigar.rbegin(), cigar.rend());

            //std::cout << cigar << std::endl;

//             printPAF((query_object->name).c_str(), query.length(), (fast_objects2.front()->name).c_str(),
//                      target.length(), k, cigar.c_str(), c, std::get<4>(region));
        }

        regions.clear();
        regions.shrink_to_fit();

        i++;
    }
}

int main(int argc, char *argv[]) {
    char optchr;
    srand(time(NULL));

//    pink::AlignmentType type = pink::global;
    int match = 2;
    int mismatch = -1;
    int gap = -2;
    std::string cigar;
//    unsigned int target_begin = 0;

    unsigned int k = 15;
    unsigned int window_length = 5;
    float f = 0.001;

    bool c = true;
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
                //type = pink::global;
                break;
            case 'S':
                //type = pink::semi_global;
                break;
            case 'L':
                //type = pink::local;
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
                c = true;
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



//        int query = rand() % fast_objects1.size();
//        int target = rand() % fast_objects1.size();
//
//        const char *q = (fast_objects1[query]->sequence).c_str();
//        const char *t = (fast_objects1[target]->sequence).c_str();
//        unsigned int q_len = (fast_objects1[query]->sequence).length();
//        unsigned int t_len = (fast_objects1[target]->sequence).length();
//
//        int cost = pink::pairwise_alignment(q, q_len, t, t_len, type, match, mismatch, gap);
//        std::cout << "\nFinal cost: " << cost << std::endl;
//
//        pink::pairwise_alignment(q, q_len, t, t_len, type, match, mismatch, gap, cigar, target_begin);
//        cigar = std::string(cigar.rbegin(), cigar.rend());
//        std::cout << "\nCigar: " << cigar << std::endl;


        //FINAL TASK!
        createQueryIndex(fast_objects1, fast_objects2, k, window_length, f, match, mismatch, gap, c);

    } else {
        printError();
    }

    return 0;
}
