#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <cctype>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <limits>
#include <ctime>
#include <utility>
#include <unordered_map>
#include <queue>
#include <deque>
#include <map>
#include <algorithm>

#include "brown_mapper.hpp"
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "brown_alignment.hpp"
#include "brown_minimizers.hpp"


const std::set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};


std::map<char, char> complement_map = {{'C', 'G'}, {'A', 'T'}, {'T', 'A'}, {'U', 'A'}, {'G', 'C'}};

// position, amount
typedef std::pair<unsigned int, unsigned int> minimizer_index_t;

// value, position, strand
typedef std::tuple<unsigned int, unsigned int, bool> triplet_t;

// query position, reference position, relative strand
typedef std::tuple<unsigned int, unsigned int, bool> minimizer_hit_t;

// starting positions, ending positions
typedef std::pair<minimizer_hit_t, minimizer_hit_t> region_t;


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
  {"threads", required_argument, NULL, 't'},
  {"band", required_argument, NULL, 'b'},
  {"cigar", no_argument, NULL, 'c'},
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
         "  -G  or  --global         global alignment (default)\n"
         "  -L  or  --local          local alignment\n"
         "  -S  or  --semi_global    semi_global alignment\n"
         "  -w  or  --window_length  <int>\n"
         "                             default: 5\n"
         "                             length of window\n"
         "  -k  or  --kmers          <int>\n"
         "                             default: 15\n"
         "                             constraints: largest supported is 16\n"
         "                             number of letters in substrings\n"
         "  -f  or  --frequency      <float>\n"
         "                             default: 0.001\n"
         "                             constraints: must be from [0, 1]\n"
         "                             number of most frequent minimizers that are not taken into account\n"
         "  -t  or  --threads        <int>\n"
         "                             default: 3\n"
         "                             number of threads\n"
         "  -b  or  --band           <int>\n"
         "                             default: 500\n"
         "                             cluster band width\n"
         "  -c  or  --cigar          include cigar string\n"
         "                             !!! will severely impact performance and size of output (quadratic memory usage and computation time) !!!\n"
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


bool hit_ordering(const minimizer_hit_t& a, const minimizer_hit_t& b) {
  if (std::get<2>(a) == std::get<2>(b)) {
    if (std::get<2>(a) == 0) {
      if (((int)std::get<0>(a) - std::get<1>(a)) == ((int)std::get<0>(b) - std::get<1>(b))) {
        return std::get<1>(a) < std::get<1>(b);
      }
      return ((int)std::get<0>(a) - std::get<1>(a)) < ((int)std::get<0>(b) - std::get<1>(b));
    } else {
      if (((int)std::get<0>(a) + std::get<1>(a)) == ((int)std::get<0>(b) + std::get<1>(b))) {
        return std::get<1>(a) < std::get<1>(b);
      }
      return ((int)std::get<0>(a) + std::get<1>(a)) > ((int)std::get<0>(b) + std::get<1>(b));
    }
  }
  return std::get<2>(a) < std::get<2>(b);
}


void prep_ref(std::vector<triplet_t>& t_minimizers, const float f) {
  std::unordered_map<unsigned int, unsigned int> ref_min_frequency;
  for (const auto& minimizer : t_minimizers) {
    ref_min_frequency[std::get<0>(minimizer)]++;
  }

  std::vector<unsigned int> occurences;
  occurences.reserve(ref_min_frequency.size());

  for (const auto& entry : ref_min_frequency) {
    occurences.push_back(entry.second);
  }

  unsigned int position = (unsigned int)((1.0f - f) * (occurences.size() - 1.0f));
  std::sort(occurences.begin(), occurences.end());

  unsigned int cutoff_freq = occurences[position] == 1 ? 2 : occurences[position];

  std::vector<triplet_t> temp;
  temp.reserve(t_minimizers.size());

  for (const auto& minimizer : t_minimizers) {
    if (ref_min_frequency[std::get<0>(minimizer)] < cutoff_freq) {
      temp.push_back(minimizer);
    }
  }

  std::swap(t_minimizers, temp);

  std::sort(t_minimizers.begin(), t_minimizers.end(),
      [] (const triplet_t& a, const triplet_t& b) {
        return (std::get<0>(a) < std::get<0>(b));
      });

  t_minimizers.shrink_to_fit();
}


std::unordered_map<unsigned int, minimizer_index_t> index_ref(
    const std::vector<triplet_t>& t_minimizers) {

  std::unordered_map<unsigned int, minimizer_index_t> ref_index;
  unsigned int pos = 0;
  unsigned int num = 0;
  unsigned int prev_min = std::get<0>(t_minimizers[0]);
  for (const auto& minimizer : t_minimizers) {
    if (prev_min != std::get<0>(minimizer)) {
      ref_index[prev_min] = std::make_pair(pos, num);
      pos += num;
      num = 1;
      prev_min = std::get<0>(minimizer); 
    } else {
      num++;
    }
  }
  ref_index[prev_min] = std::make_pair(pos, num);

  return ref_index;
}


std::vector<minimizer_hit_t> find_minimizer_hits(
    const std::unordered_map<unsigned int, minimizer_index_t>& ref_index,
    const std::vector<triplet_t>& t_minimizers,
    const std::vector<triplet_t>& q_minimizers) {

  std::vector<minimizer_hit_t> hits;

  for (const auto& minimizer : q_minimizers) {
    auto found = ref_index.find(std::get<0>(minimizer));
    if (found != ref_index.end()) {
      for (unsigned int i = 0; i < found->second.second; i++) {
        if (std::get<2>(minimizer) == std::get<2>(t_minimizers[found->second.first + i])) {
          hits.emplace_back(
            std::get<1>(minimizer),
            std::get<1>(t_minimizers[found->second.first + i]),
            0);  
        } else {
          hits.emplace_back(
            std::get<1>(minimizer),
            std::get<1>(t_minimizers[found->second.first + i]),
            1);
        }
        
      }
    }
  }

  return hits;
}


int band(const minimizer_hit_t& h2, const minimizer_hit_t& h1) {
  if (std::get<2>(h2) == 0) {
    return ((int)std::get<0>(h2) - (int)std::get<1>(h2)) - ((int)std::get<0>(h1) - (int)std::get<1>(h1));
  }
  return ((int)std::get<0>(h1) + (int)std::get<1>(h1)) - ((int)std::get<0>(h2) + (int)std::get<1>(h2));
}


bool region_overlap(const region_t& a, const region_t&b, int bw) {
  if (std::get<2>(a.first) == std::get<2>(b.first)) {
    if (std::get<2>(a.first) == 0) {
      return std::abs(band(a.second, b.first)) < bw;
    } else {
      return std::abs(band(std::make_tuple(std::get<0>(a.second), std::get<1>(a.first), 1), std::make_tuple(std::get<0>(b.first), std::get<1>(b.second), 1))) < bw;
    }
  }
  return false;
}


void merge_regions(std::vector<std::pair<region_t, unsigned int>>& regions, int bw) {
  if (regions.size() < 1) {
    return;
  }
  std::sort(regions.begin(), regions.end(),
      [] (const std::pair<region_t, unsigned int>& a,
          const std::pair<region_t, unsigned int>& b) {
        return std::get<1>(a.first.first) > std::get<1>(b.first.first);
      });
  unsigned int index = 0;
  for (unsigned int i = 0; i < regions.size(); ++i) {
    if (index > 0 && region_overlap(regions[index - 1].first, regions[i].first, bw)) {
      while (index > 0 && region_overlap(regions[index - 1].first, regions[i].first, bw)) {
        regions[index - 1].first.second = std::make_tuple(std::max(std::get<0>(regions[index - 1].first.second), std::get<0>(regions[i].first.second)),
                                                    std::max(std::get<1>(regions[index - 1].first.second), std::get<1>(regions[i].first.second)),
                                                    std::get<2>(regions[i].first.second));
        regions[index - 1].first.first = std::make_tuple(std::min(std::get<0>(regions[index - 1].first.first), std::get<0>(regions[i].first.first)),
                                                   std::min(std::get<1>(regions[index - 1].first.first), std::get<1>(regions[i].first.first)),
                                                   std::get<2>(regions[i].first.second));
        regions[index - 1].second += std::max(0, (int)std::get<1>(regions[i].first.second) - (int)std::get<1>(regions[index - 1].first.second));
        index--;
      }
    } else {
      regions[index] = regions[i];
    }
    index++;
  }
  regions.erase(regions.begin() + index, regions.end());
}


std::pair<region_t, unsigned int> find_region(std::vector<minimizer_hit_t>& hits,
    unsigned int b, unsigned int e, unsigned int k) {

  std::pair<region_t, unsigned int> region;
  region = std::make_pair(std::make_pair(std::make_tuple(0, 0, 0), std::make_tuple(0, 0, 0)), 0);
  if (hits.empty()) {
    return region;
  }

  std::sort(hits.begin() + b, hits.begin() + e + 1, 
      [] (const minimizer_hit_t& a, const minimizer_hit_t& b) noexcept {
        if (std::get<2>(a) == 0) {
          if (std::get<0>(a) == std::get<0>(b)) {
            return std::get<1>(a) < std::get<1>(b);
          }
          return std::get<0>(a) < std::get<0>(b);
        } else {
          if (std::get<0>(a) == std::get<0>(b)) {
            return std::get<1>(a) < std::get<1>(b);
          }
          return std::get<0>(a) > std::get<0>(b);
        }
      });

  std::vector<region_t> lis(e - b + 1,
      std::make_pair(std::make_tuple(0, 0, 0), std::make_tuple(0, 0, 0)));

  std::vector<unsigned int> matches(e - b + 1, 0);

  unsigned int len = 1;
  lis[0] = std::make_pair(hits[b], hits[b]);
  matches[0] = k;

  for (unsigned int i = b + 1; i < e + 1; ++i) {
    if (std::get<1>(hits[i]) > std::get<1>(lis[len - 1].second)) {
      unsigned int diff = std::get<0>(hits[i]) - std::get<0>(lis[len - 1].second);
      diff = std::min(diff, std::get<1>(hits[i]) - std::get<1>(lis[len - 1].second));
      diff = std::min(diff, k);
      matches[len] = matches[len - 1] + diff;
      lis[len] = lis[len - 1];
      lis[len].second = hits[i];
      len++;
    } else {
      auto pair_it = std::upper_bound(lis.begin(), lis.begin() + len, hits[i],
          [] (const minimizer_hit_t& a,
              const region_t& b) {
            return (std::get<1>(a) < std::get<1>(b.second));
          });
      if (pair_it == lis.begin()) {
        *pair_it = std::make_pair(hits[i], hits[i]);
      } else {
        *pair_it = std::make_pair((pair_it - 1)->first, hits[i]);
        unsigned int diff = std::get<0>(hits[i]) - std::get<0>((pair_it - 1)->second);
        diff = std::min(diff, std::get<1>(hits[i]) - std::get<1>((pair_it - 1)->second));
        diff = std::min(diff, k);
        matches[pair_it - lis.begin()] = matches[pair_it - 1 - lis.begin()] + diff;
      }
    }
  }
  if (len < 1) {
    return region;
  }
  region_t max_reg = lis[len - 1];
  if (std::get<2>(max_reg.first) == 1) {
    unsigned int temp = std::get<0>(max_reg.first);
    std::get<0>(max_reg.first) = std::get<0>(max_reg.second);
    std::get<0>(max_reg.second) = temp;
  }
  region = std::make_pair(max_reg, matches[len - 1]);
  return region; 
}


unsigned int read_cigar(const std::string& cigar, unsigned int& size) {
  size = 0;
  char* next;
  unsigned int matches = 0;
  unsigned int pos = 0;
  while (pos < cigar.size()) {
    unsigned int num = strtol(cigar.c_str() + pos, &next, 10);
    pos = next - cigar.c_str();
    if (cigar[pos] == '=') {
      matches += num;
    } else if (cigar[pos] == 'X' || cigar[pos] == 'I' || cigar[pos] == 'D') {
      size += num;
    }
    pos++;
  }
  size += matches;
  return matches;
}


std::string reverse_complement(const std::string& original, unsigned int pos, unsigned int length) {
  std::string rc(original.begin() + pos, original.begin() + pos + length);
  unsigned int j = pos + length - 1;
  for (unsigned int i = 0; i < length; ++i) {
    rc[i] = complement_map[original[j--]];
  }
  return rc;
}

// Sort by strand, diff, tpos
// Cluster
// For each cluster:
//   Sort by position
//   LIS
//   Store regions
//   Merge overlaping regions | ones in epsilon range
std::string map_paf(const std::vector<triplet_t>& t_minimizers,
    const std::unordered_map<unsigned int, minimizer_index_t>& ref_index,
    const std::vector<std::unique_ptr<FastAQ>>& fastaq_objects1,
    const std::vector<std::unique_ptr<FastAQ>>& fastaq_objects2,
    unsigned int k, unsigned int window_length, brown::AlignmentType alignment,
    int match, int mismatch, int gap,
    unsigned int t_begin, unsigned int t_end, int bw, bool c) {

  std::string paf;
  for (unsigned int fobj = t_begin; fobj < t_end; fobj++) {
    std::vector<triplet_t> q_minimizers = brown::minimizers(
      fastaq_objects1[fobj]->sequence.c_str(), fastaq_objects1[fobj]->sequence.size(), k, window_length);

    std::vector<minimizer_hit_t> hits = find_minimizer_hits(ref_index, t_minimizers, q_minimizers);

    std::sort(hits.begin(), hits.end(), hit_ordering);

    std::vector<std::pair<region_t, unsigned int>> regions;

    unsigned int b = 0;
    for (unsigned int e = 0; e < hits.size(); ++e) {
      if (e == hits.size() - 1 ||
          std::get<2>(hits[e + 1]) != std::get<2>(hits[e]) ||
          std::abs(band(hits[e + 1], hits[e])) >= bw) {
        std::pair<region_t, unsigned int> temp_region = find_region(hits, b, e, k);
        if (temp_region.second > 4 * k) {
          regions.push_back(temp_region);
        }
        b = e + 1;
      }
    }

    if (regions.empty()) {
      continue;
    }

    merge_regions(regions, bw);

    std::sort(regions.begin(), regions.end(),
        [] (const std::pair<region_t, unsigned int>& a,
            const std::pair<region_t, unsigned int>& b) {
      return a.second > b.second;
    });

    for (const auto& reg : regions) {
      region_t max_region = reg.first;
      unsigned int matches = reg.second;

      std::string rel_strand;
      unsigned int total;
      unsigned int q_b, q_e, t_b, t_e;
      std::string cigar;
      unsigned int target_begin;

      t_b = std::get<1>(max_region.first);
      t_e = std::get<1>(max_region.second) + k;
      q_b = std::get<0>(max_region.first);
      q_e = std::get<0>(max_region.second) + k;
      total = (t_e - t_b);
      
      if (!std::get<2>(max_region.first)) {
        rel_strand = "+";
        if (c) {
          brown::pairwise_alignment(fastaq_objects1[fobj]->sequence.c_str() + q_b,
              q_e - q_b,
              fastaq_objects2[0]->sequence.c_str() + t_b,
              t_e - t_b,
              alignment, match, mismatch, gap, cigar, target_begin);
          matches = read_cigar(cigar, total);
        }
      } else {
        rel_strand = "-";
        if (c) {
          std::string rc = reverse_complement(fastaq_objects1[fobj]->sequence, q_b,
              q_e - q_b);
          brown::pairwise_alignment(rc.c_str(), rc.size(),
              fastaq_objects2[0]->sequence.c_str() + t_b,
              t_e - t_b,
              alignment, match, mismatch, gap, cigar, target_begin);
          matches = read_cigar(cigar, total);
        }
      }
      paf += fastaq_objects1[fobj]->name + "\t" +
             std::to_string(fastaq_objects1[fobj]->sequence.size()) + "\t" +
             std::to_string(q_b) + "\t" +
             std::to_string(q_e) + "\t" +
             rel_strand + "\t" +
             fastaq_objects2[0]->name + "\t" +
             std::to_string(fastaq_objects2[0]->sequence.size()) + "\t" +
             std::to_string(t_b) + "\t" +
             std::to_string(t_e) + "\t" +
             std::to_string(matches) + "\t" +
             std::to_string(total) + "\t" +
             std::to_string(255);
      if (c) {
        paf += "\tcg:Z:" + cigar;
      }
      paf += "\n";
    }
  }
  return paf;
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
  float f = 0.001f;

  unsigned int t = 3;
  int bw = 500;
  bool c = false;

  while ((optchr = getopt_long(argc, argv, "hvm:g:M:GLSk:w:f:t:b:c", long_options, NULL)) != -1) {
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
        f = atof(optarg);
        if (f < 0.0 || f > 1.0) {
          fprintf(stderr, "[mapper] error: f must be from [0, 1]!\n"); 
          exit(1); 
        }
        break;
      }
      case 't': {
        t = atoi(optarg);
        if (t < 1) {
          fprintf(stderr, "[mapper] error: t must be positive!\n"); 
          exit(1); 
        }
        break;
      }
      case 'b': {
        bw = atoi(optarg);
        if (bw < 1) {
          fprintf(stderr, "[mapper] error: b must be positive!\n"); 
          exit(1);
        }
        break;
      }
      case 'c': {
        c = true;
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

  // Load files
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

  fprintf(stderr, " Done!\n\nParsing files...");

  // Parse files

  std::vector<std::unique_ptr<FastAQ>> fastaq_objects1;
  std::vector<std::unique_ptr<FastAQ>> fastaq_objects2;

  FastAQ::parse(fastaq_objects1, file1, file1_format);
  FastAQ::parse(fastaq_objects2, file2, file2_format);

  fprintf(stderr, " Done!\n");

  // Print file stats
  // fprintf(stderr, "\n");  
  // FastAQ::print_statistics(fastaq_objects1, file1);
  // FastAQ::print_statistics(fastaq_objects2, file2);

  // Index reference and map
  fprintf(stderr, "\nMapping process started with parameters:\n"
                  "  k             = %u\n"
                  "  window length = %u\n"
                  "  top f freq    = %g\n"
                  "  band width    = %u\n"
                  "  threads       = %u\n",
                  k, window_length, f, bw, t);
  if (c) {
    fprintf(stderr, " Alignment\n"
                    "  type     = %s\n"
                    "  match    = %d\n"
                    "  mismatch = %d\n"
                    "  gap      = %d\n",
                    alignmentType.c_str(), match, mismatch, gap);
  }
  fprintf(stderr, " Indexing reference...");

  // Parallelized version of reference minimizers search
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_ref = thread_pool::createThreadPool(t);

  std::vector<std::future<std::vector<triplet_t>>> thread_futures_ref;

  for (unsigned int tasks = 0; tasks < t - 1; ++tasks) {
    thread_futures_ref.emplace_back(thread_pool_ref->submit_task(brown::minimizers,
        fastaq_objects2[0]->sequence.c_str() + tasks * fastaq_objects2[0]->sequence.size() / t,
        fastaq_objects2[0]->sequence.size() / t - 1 + window_length + k - 1,
        k, window_length));
  }
  thread_futures_ref.emplace_back(thread_pool_ref->submit_task(brown::minimizers,
        fastaq_objects2[0]->sequence.c_str() + (t - 1) * fastaq_objects2[0]->sequence.size() / t,
        fastaq_objects2[0]->sequence.size() - (t - 1) * fastaq_objects2[0]->sequence.size() / t,
        k, window_length));

  std::vector<triplet_t> t_minimizers;
  for (unsigned int i = 0; i < t; ++i) {
    thread_futures_ref[i].wait();
    unsigned int offset = i * fastaq_objects2[0]->sequence.size() / t;
    for (auto& el : thread_futures_ref[i].get()) {
      std::get<1>(el) += offset;
      t_minimizers.push_back(el);
    }
  }

  // Regular version of reference minimizers search
  // std::vector<triplet_t> t_minimizers = brown::minimizers(
  //     fastaq_objects2[0]->sequence.c_str(), fastaq_objects2[0]->sequence.size(), k, window_length);

  prep_ref(t_minimizers, f);

  std::unordered_map<unsigned int, minimizer_index_t> ref_index = index_ref(t_minimizers);

  fprintf(stderr, " Done!\n Mapping...");

  std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(t);

  std::vector<std::future<std::string>> thread_futures;

  for (unsigned int tasks = 0; tasks < t - 1; ++tasks) {
    thread_futures.emplace_back(thread_pool->submit_task(map_paf, std::ref(t_minimizers), std::ref(ref_index),
      std::ref(fastaq_objects1), std::ref(fastaq_objects2),
      k, window_length, alignment, match, mismatch, gap, tasks * fastaq_objects1.size() / t, (tasks + 1) * fastaq_objects1.size() / t, bw, c));  
  }
  thread_futures.emplace_back(thread_pool->submit_task(map_paf, std::ref(t_minimizers), std::ref(ref_index),
      std::ref(fastaq_objects1), std::ref(fastaq_objects2),
      k, window_length, alignment, match, mismatch, gap, (t - 1) * fastaq_objects1.size() / t, fastaq_objects1.size(), bw, c));
  
  for (auto& it : thread_futures) {
    it.wait();
    std::cout << it.get();
  }

  fprintf(stderr, " Done!\nFinished mapping process!\n");

  return 0;
}
