#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <limits>
#include <ctime>
#include <unordered_map>
#include <algorithm>

#include "brown_mapper.hpp"
#include "bioparser/bioparser.hpp"
#include "brown_alignment.hpp"
#include "brown_minimizers.hpp"

const std::set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};

static struct option long_options[] = {
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {"match", required_argument, NULL, 'M'},
  {"mismatch", required_argument, NULL, 'm'},
  {"gap", required_argument, NULL, 'g'},
  {"local", no_argument, NULL, 'L'},
  {"global", no_argument, NULL, 'G'},
  {"semi_global", no_argument, NULL, 'S'},
  {"window_length", required_argument, NULL, 'w'},
  {"kmers", required_argument, NULL, 'k'},
  {"frequency", required_argument, NULL, 'f'},
  {NULL, no_argument, NULL, 0}
};

class FastAQ {
public:
  std::string name;
  std::string sequence;
  std::string quality;

  FastAQ(
    const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length) : FastAQ(name, name_length, sequence, sequence_length, "", 0){}

  FastAQ(
    const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length,
    const char* quality, uint32_t quality_length) {
      this->name = {name, name_length};
      this->sequence = {sequence, sequence_length};
      this->quality = {quality, quality_length};
  }

  static void parse(
    std::vector<std::unique_ptr<FastAQ>> &fastaq_objects,
    const std::string file, const int file_format) {
      if (file_format == 1) {
        auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FastAQ>(file);
        fasta_parser->parse_objects(fastaq_objects, -1);
      } else {
        auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FastAQ>(file);
        fastq_parser->parse_objects(fastaq_objects, -1);
      }
    }

  static void print_statistics(
    const std::vector<std::unique_ptr<FastAQ>> &fastaq_objects,
    const std::string file) {
      int num = fastaq_objects.size();
      double average = 0;
      uint32_t max = 0;
      uint32_t min = std::numeric_limits<int>::max();
      for (int i = 0; i < num; i++) {
        average += fastaq_objects[i]->sequence.size();
        if (fastaq_objects[i]->sequence.size() > max) {
          max = fastaq_objects[i]->sequence.size();
        }
        if (fastaq_objects[i]->sequence.size() < min) {
          min = fastaq_objects[i]->sequence.size();
        }
      }
      average /= num;
      fprintf(stderr, "Stats for: %s\n"
                      "  Number of sequences: %d\n"
                      "  Average length:      %g\n"
                      "  Maximum length:      %d\n"
                      "  Minimum length:      %d\n",
                      file.c_str(), num, average, max, min);
  }
};

void help(void) {
  printf("brown_mapper - for mapping long erroneous fragments from third\n"
         "               generation of sequencing technologies to a reference\n"
         "               genome.\n\n"

         "Usage: brown_mapper [OPTIONS] [file1 file2]   start mapper\n"
         "file1 - FASTA/FASTQ file containing a set of fragments\n"
         "file2 - FASTA file containing reference genome\n\n"

         "Supported file extensions: .fasta\n"
         "                           .fa\n"
         "                           .fastq\n"
         "                           .fq\n"
         "                           .fasta.gz\n"
         "                           .fa.gz\n"
         "                           .fastq.gz\n"
         "                           .fq.gz\n\n"
         "OPTIONS:\n"
         "  -h  or  --help           print help (displayed now) and exit\n"
         "  -v  or  --version        print version info and exit\n"
         "  -M  or  --match          <int>\n"
         "                             default: 4\n"
         "                             match number\n"
         "  -m  or  --mismatch       <int>\n"
         "                             default: -1\n"
         "                             mismatch number\n"
         "  -g  or  --gap            <int>\n"
         "                             default: -2\n"
         "                             gap number\n"
         "  -G  or  --global         global alignment\n"
         "  -L  or  --local          local alignment\n"
         "  -S  or  --semi_global    semi_global alignment\n"
         "  -w  or  --window_length  <int>\n"
         "                             default: 5\n"
         "                             length of window\n"
         "  -k  or  --kmers          <int>\n"
         "                             default: 15\n"
         "                             number of letters in substrings\n"
         "  -f  or  --frequency      <int>\n"
         "                             default: 0.001\n"
         "                             number of frequent minimizers that are not taken in account\n"
  );
}

void version(void) {
  printf("brown_mapper %d.%d\n",
    brown_mapper_VERSION_MAJOR,
    brown_mapper_VERSION_MINOR
  );
}

bool contains_extension(const std::string file, const std::set<std::string> &extensions) {
  for (const auto& it: extensions) {
    if (file.size() > it.size()) {
      if (file.compare(file.size()-it.size(), std::string::npos, it) == 0) {
        return true;
      }
    }
  }
  return false;
}

int main (int argc, char **argv) {
  int optchr;

  int match = 4;
  int mismatch = -1;
  int gap = -2;
  brown::AlignmentType alignment = brown::AlignmentType::global;
  std::string alignmentType = "global";

  int window_length = 5;
  int k = 15;
  float f = 0.001;

  while ((optchr = getopt_long(argc, argv, "hvm:g:M:GLSk:w:f:", long_options, NULL)) != -1) {
    switch (optchr) {
      case 'h': {
        help();
        exit(0);
      }
      case 'v': {
        version();
        exit(0);
      }
      case 'M': {
        match = atoi(optarg);
        break;
      }
      case 'm': {
        mismatch = atoi(optarg);
        break;
      }
      case 'g': {
        gap = atoi(optarg);
        break;
      }
      case 'G': {
        alignment = brown::AlignmentType::global;
        alignmentType = "global";
        break;
      }
      case 'L': {
        alignment = brown::AlignmentType::local;
        alignmentType = "local";
        break;
      }
      case 'S': {
        alignment = brown::AlignmentType::semi_global;
        alignmentType = "semi_global";
        break;
      }
      case 'w': {
        window_length = atoi(optarg);
        break;
      }
      case 'k': {
        k = atoi(optarg);
        break;
      }
      case 'f': {
        f = atoi(optarg);
        break;
      }
      default: {
        fprintf(stderr, "[mapper] error: Unknown option. Type %s --help for usage.\n", argv[0]);
        exit(1);
      }
    }
  }

  if (argc - optind != 2) {
    fprintf(stderr, "[mapper] error: Expected 2 mapping arguments! Use --help for usage.\n");
    exit(1);
  }

  fprintf(stderr, "\nLoading files...");

  std::string file1 (argv[optind]);
  std::string file2 (argv[optind+1]);

  int file1_format = 0;
  int file2_format = 0;

  if (contains_extension(file1, fasta_formats)) {
    file1_format = 1;
  } else if (contains_extension(file1, fastq_formats)) {
    file1_format = 2;
  }

  file2_format = contains_extension(file2, fasta_formats);

  if ( !(file1_format && file2_format) ) {
    fprintf(stderr, "[mapper] error: Unsupported format(s)! Check --help for supported file formats.\n");
    exit(1);
  }

  fprintf(stderr, " Done!\n");
  fprintf(stderr, "\nParsing files...");

  std::vector<std::unique_ptr<FastAQ>> fastaq_objects1;
  std::vector<std::unique_ptr<FastAQ>> fastaq_objects2;

  FastAQ::parse(fastaq_objects1, file1, file1_format);
  FastAQ::parse(fastaq_objects2, file2, file2_format);

  fprintf(stderr, " Done!\n\n");



  // File stats

  FastAQ::print_statistics(fastaq_objects1, file1);
  FastAQ::print_statistics(fastaq_objects2, file2);



  // Minimizers

  fprintf(stderr, "\nFinding minimizers with parameters:\n");
  fprintf(stderr, "  k = %d\n  Window length = %d\n  f= %f\n", k, window_length, f);

  char prog[] = {'|', '/', '-', '\\'};

  int counter = 0;
  int charcount = 0;
  int percentage = 0.01 * fastaq_objects1.size();
  std::unordered_map<unsigned int, unsigned int> frequency_map;
  for(unsigned int i = 0; i < fastaq_objects1.size(); i++){

  	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = brown::minimizers(fastaq_objects1[i]->sequence.c_str(),
  																				fastaq_objects1[i]->sequence.size(),k , window_length);
  	for(auto const& minimizer: minimizers){
  		frequency_map[std::get<0>(minimizer)]++;
  	}
  	counter++;
    if(counter % percentage == 0){
      charcount++;
      fprintf(stderr, "\r[ %c %3.0f%% ]", prog[charcount%4], (float)counter/fastaq_objects1.size() * 100);
    }
  }

  fprintf(stderr, "\rMinimizers found!\n");

  fprintf(stderr, "Number of distinct minimizers: %lu\n", frequency_map.size());

  std::vector<unsigned int> occurences;
  occurences.reserve(frequency_map.size());

  for(auto const& entry : frequency_map){
  	occurences.push_back(entry.second);
  }

  int position = (1-f) * occurences.size();

  std::sort(occurences.begin(), occurences.end());

  fprintf(stderr, "Number of occurences of the most frequent minimizer: %d\n", occurences[position]);

  fprintf(stderr, " Done!\n");



  // Alignment

  // fprintf(stderr, "\nStarting alignment with parameters:\n");
  // fprintf(stderr, "  Match = %d \n  Mismatch = %d\n  Gap/Indel = %d\n", match, mismatch, gap);
  // std::cerr << "  Alignment type = " << alignmentType <<std::endl;
  // fprintf(stderr, "Aligning...");

  // srand(time(NULL));
  // int i1 = rand() % fastaq_objects1.size();
  // int i2 = rand() % fastaq_objects1.size();

  // std::string cigar;
  // unsigned int target_begin = 0;
  // int value = brown::pairwise_alignment(fastaq_objects1[i1]->sequence.c_str(), fastaq_objects1[i1]->sequence.size(),
  //                                       fastaq_objects1[i2]->sequence.c_str(), fastaq_objects1[i2]->sequence.size(),
  //                                       alignment, match, mismatch, gap, cigar, target_begin);
  // std::cout << "Alignment score: " << value << std::endl;
  // std::cout << "Pos: " << target_begin << std::endl;
  // std::cout << "CIGAR: " << cigar << std::endl;

  // fprintf(stderr, " Done!\n");
  


  // Small test examples

  // std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = brown::minimizers("CGAC",4, 3 , 1);
  // for(auto const& minimizer : minimizers){
  // 	   std::cout << "minimizer: " << std::get<0>(minimizer) << std::endl;
  // 	   std::cout << "pozicija: " << std::get<1>(minimizer) << std::endl;
  // 	   std::cout << "bool: " << std::get<2>(minimizer) << std::endl;
  // }

  // std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = brown::minimizers("CGACT", 5, 3, 1);
  // for(auto const& minimizer : minimizers){
  // 	   std::cout << "minimizer: " << std::get<0>(minimizer) << std::endl;
  // 	   std::cout << "pozicija: " << std::get<1>(minimizer) << std::endl;
  // 	   std::cout << "bool: " << std::get<2>(minimizer) << std::endl;
  // }

  // std::string q = {"ATTGGAA"};
  // std::string t = {"CCCACTTTTTGG"};
  // std::string cigar;
  // unsigned int target_begin = 0;
  // int value = brown::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), brown::AlignmentType::semi_global, match, mismatch, gap, cigar, target_begin);
  // std::cout << value << std::endl;
  // std::cout << target_begin << std::endl;
  // std::cout << cigar << std::endl;

  // std::string q = {"TTCCGCCAA"};
  // std::string t = {"AACCCCTT"};
  // std::string cigar;
  // unsigned int target_begin = 0;
  // int value = brown::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), brown::AlignmentType::local, match, mismatch, gap, cigar, target_begin);
  // std::cout << value << std::endl;
  // std::cout << target_begin << std::endl;
  // std::cout << cigar << std::endl;

  // std::string q = {"TGCATAT"};
  // std::string t = {"ATCCGAT"};
  // std::string cigar;
  // unsigned int target_begin = 0;
  // int value = brown::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), brown::AlignmentType::global, match, mismatch, gap, cigar, target_begin);
  // std::cout << value << std::endl;
  // std::cout << target_begin << std::endl;
  // std::cout << cigar << std::endl;

  fprintf(stderr, "\nFinished!\n");

  return 0;
}
