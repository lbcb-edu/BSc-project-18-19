#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <limits>

#include "brown_mapper.hpp"
#include "bioparser/bioparser.hpp"

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
      this->name = name;
      this->sequence = sequence;
      this->quality = quality;
  }

  static void parse(
    std::vector<std::unique_ptr<FastAQ>> &fastaq_objects,
    const std::string &file, const int &file_format) {
      if (file_format == 1) {
        auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FastAQ>(file);
        fasta_parser->parse_objects(fastaq_objects, -1);
      } else {
        auto fastq_parser = bioparser::createParser<bioparser::FastqParser, FastAQ>(file);
        uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
        while (true) {
          auto status = fastq_parser->parse_objects(fastaq_objects, size_in_bytes);
          if (status == false) {
            break;
          }
        }
      }
      FastAQ::print_statistics(fastaq_objects, file);
    }

    static void print_statistics(
      const std::vector<std::unique_ptr<FastAQ>> &fastaq_objects,
      const std::string &file) {
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
         "               genome, which has various use cases in bioinformatics.\n\n"

         "Usage: brown_mapper [OPTIONS] [file1 file2]   start mapping files\n"
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
         "  -h  or  --help         print help (displayed now) and exit\n"
         "  -v  or  --version      print version info and exit\n"
  );
}

void version(void) {
  printf("brown_mapper %d.%d\n",
    brown_mapper_VERSION_MAJOR,
    brown_mapper_VERSION_MINOR
  );
}

bool contains_extension(const std::string &file, const std::set<std::string> &extensions) {
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
  static struct option long_options[] = {
    {"help",    no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    {NULL,      no_argument, NULL,  0 }
  };

  while ((optchr = getopt_long(argc, argv, "hv", long_options, NULL)) != -1) {
    switch (optchr) {
      case 'h': {
        help();
        exit(0);
      }
      case 'v': {
        version();
        exit(0);
      }
      default: {
        printf("Unknown option. Type %s --help for usage.\n", argv[0]);
        exit(1);
      }
    }
  }

  if (argc - optind != 2) {
    printf("Expected 2 mapping arguments! Use --help for usage.\n");
    exit(1);
  }

  std::set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
  std::set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};

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
    printf("Unsupported format(s)! Check --help for supported file formats.\n");
    exit(1);
  }

  std::vector<std::unique_ptr<FastAQ>> fastaq_objects1;
  std::vector<std::unique_ptr<FastAQ>> fastaq_objects2;

  FastAQ::parse(fastaq_objects1, file1, file1_format);
  FastAQ::parse(fastaq_objects2, file2, file2_format);

  return 0;
}
