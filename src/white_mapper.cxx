#include <iostream>
#include <memory>
#include <vector>
#include <getopt.h>
#include <stdlib.h>
#include <set>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <future>
#include <sstream>

#include "white_minimizers.h"
#include "white_alignment.h"
#include "white_mapper.h"
#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"

#define MATCH 2;
#define MISMATCH -1;
#define GAP -2;
#define K_MER_LENGTH 15;
#define WINDOW_SIZE 5;
#define PERCENTAGE 0.001;
#define NUMBER_OF_THREADS 5;

class SequenceFormat
{
	 public:
                std::string name;
                std::string sequence;
		std::string quality;

		friend bioparser::FastaParser<SequenceFormat>;
		friend bioparser::FastqParser<SequenceFormat>;
		friend void statistics (std::vector<std::unique_ptr<SequenceFormat>>&);

	private:
		SequenceFormat (
                        const char* name, uint32_t name_length,
                        const char* sequence, uint32_t sequence_length,
			const char* quality, uint32_t quality_length
			) : name (name, name_length), sequence (sequence, sequence_length), quality (quality, quality_length)
                        {
			}

		SequenceFormat (
                        const char* name, uint32_t name_length,
                        const char* sequence, uint32_t sequence_length
                        ) : name (name, name_length), sequence (sequence, sequence_length), quality ("")
                        {
                        }

};

const std::set<std::string> fasta_formats = {".fasta", ".fa", ".fasta.gz", ".fa.gz"};
const std::set<std::string> fastq_formats = {".fastq", ".fq", ".fastq.gz", ".fq.gz"};


const struct option long_options[] = {
{"help",	no_argument,	0,	'h' },
{"version",	no_argument,	0,	'v' },
{0,		0,		0,	 0  },
};

void statistics (std::vector<std::unique_ptr<SequenceFormat>> &v, std::string file)
{
        uint32_t min_length = v.at(0) -> sequence.size();
        uint32_t max_length = min_length;
        uint32_t total_length = 0;

        for (auto &ptr : v)
        {
                if (ptr -> sequence.size() < min_length)
                        min_length = ptr -> sequence.size();

                if (ptr -> sequence.size() > max_length)
                        max_length = ptr -> sequence.size();

                total_length += ptr -> sequence.size();
        }

	printf( "\nStatistics for file: %s\n"
		"	Number of sequences: %ld\n"
		"	Minimum length: %d\n"
		"	Maximum length: %d\n"
		"	Average length: %.2f\n\n" ,
		file.c_str(), v.size(), min_length, max_length, (float)total_length/v.size() );

}

void help() {
	std::cout << 		"--Help:"
				"\n	This program is a mapper that requires two files. First file should contain a set of	\n"
				"	fragments in FASTA or FASTQ format and the second one a corresponding reference genome	\n"
                                "	in FASTA format. The mapper parses the files and stores them in	memory while also	\n"
				"	providing statistics of given files:							\n"
                                "	-> number of sequences									\n"
                                "	-> average length									\n"
                                "	-> minimal and maximal length.								\n\n"

				"	Afterwards, two random sequences from the first file are selected and aligned using	\n"
				"	global, semi-global and local alignment. The resulting cigar string is printed.		\n"
				"	Default for match, mismatch ang gap values is 2, -1, -2 respectively. The default	\n"
				"	values can be changed using appropriate options.					\n\n"

				"	Seed and extend:									\n"
				"	In order to alleviate the alignment proces we use seed and extend approach. Among all	\n"
				"	k-mers we choose a set of minimizers which will be used for fast detection of similar	\n"
				"	regions between two sequences prior the exact alignment. The mapper will find minimizers\n"
				"	for every sequence in the first file and print the number of distinct minimizers and the\n"
				"	number of occurences of the most frequent minimizer when the top f frequent minimizers 	\n"
				"	are not taken in account.								\n\n"

				"	The proces of finding minimizers uses three variables:					\n"
				"	-> k-mer length k									\n"
				"	-> window size w									\n"
				"	-> percentage of top minimizers to disregard f						\n\n"

				"	Their default values are 15, 5 and 0.001 respecively.					\n\n"

                                "	File extensions accepted:								\n"
                                "	-> .fasta             -> .fastq								\n"
                                "	-> .fa                -> .fq								\n"
                                "	-> .fasta.gz          -> .fastq.gz							\n"
                                "	-> .fa.gz             -> .fq.gz								\n\n"

				"	Usage: white_mapper [OPTIONS] sequences_file reference_genome_file match mismatch gap	\n\n"

                                "	There are five options:									\n"
                                "	-> \"-h\" or \"--help\" for help							\n"
                                "	-> \"-v\" or \"--version\" for displaying the current version.				\n"
				"	-> \"-m ARG\" sets match value to ARG							\n"
				"       -> \"-s ARG\" sets mismatch value to ARG                                                \n"
				"       -> \"-g ARG\" sets gap value to ARG                                                     \n"
				"       -> \"-k ARG\" sets length k of k-mers to ARG						\n"
				"       -> \"-w ARG\" sets window size to ARG							\n"
				"	-> \"-f ARG\" sets f to ARG								\n\n";
}

void version() {
	std::cout << "\n--Version:";
	std::cout << " " <<white_mapper_VERSION_MAJOR <<"." <<white_mapper_VERSION_MINOR <<"\n\n";
}

bool check_format (std::string file, std::set<std::string> list_of_formats) {
	for (auto &format : list_of_formats)
		if(file.size() > format.size())
			if (file.compare (file.size() - format.size(), format.size(), format) == 0)
				return true;
	return false;
}

template < template<class> class T>
std::vector<std::unique_ptr<SequenceFormat>> parse_file (std::string file_path)
{
	std::vector<std::unique_ptr<SequenceFormat>> objects;
	auto parser = bioparser::createParser<T, SequenceFormat> (file_path);
	parser -> parse_objects (objects, -1);

	statistics (objects, file_path);

	return objects;
}

bool comp(std::vector< std::tuple <unsigned int, unsigned int, bool>> a, std::vector <std::tuple <unsigned int, unsigned int, bool>> b) {
	return a.size() < b.size();
}

std::vector<std::tuple<unsigned int, unsigned int, bool>> reference_genome_minimizers (std::unique_ptr<SequenceFormat> &ref, unsigned int k, unsigned int window_size, float f) {
	std::unordered_map <unsigned int, unsigned int > minimizer_position_in_vector;
	std::vector <std::tuple <unsigned int, unsigned int, bool>> minimizers;
	std::vector <std::vector <std::tuple<unsigned int, unsigned int, bool>>> occurrences;

	minimizers = white::minimizers(ref -> sequence.c_str(), ref -> sequence.size(), k, window_size);

	unsigned int i = 0;

	std::vector <std::tuple <unsigned int, unsigned int, bool>> temp;

	for (auto minimizer_tuple : minimizers) {
		if (minimizer_position_in_vector[std::get<0>(minimizer_tuple)] == 0){
			i++;
			temp.emplace_back(minimizer_tuple);
			minimizer_position_in_vector[std::get<0>(minimizer_tuple)] = i;
			occurrences.emplace_back(temp);
			temp.clear();
		}
		else
			occurrences[minimizer_position_in_vector[std::get<0>(minimizer_tuple)] - 1].emplace_back(minimizer_tuple);
	}

	unsigned int number_of_minimizers_to_keep = (unsigned int) std::round(occurrences.size() * (1 - f));

	std::sort(occurrences.begin(), occurrences.end(), comp);

	std::vector <std::tuple <unsigned int, unsigned int, bool>> result;

	for (unsigned int i = 0; i < number_of_minimizers_to_keep; i++)
		for (auto tuple : occurrences[i])
			result.emplace_back(tuple);

	return result;

}

bool compare_minimizers(std::tuple <unsigned int, unsigned int, bool> a, std::tuple <unsigned int, unsigned int, bool> b) {
	return std::get<0>(a) < std::get<0>(b);
}

int my_binary_search(std::vector <std::tuple <unsigned int, unsigned int, bool>> &array, unsigned int number_to_find) {
	int min = 0;
	int max = (array.size() - 1);
	int guess;
	int result = -1;

	while (min <= max)
	{
		guess = (int)(((max + min) / 2) + 0.5);

		if (number_to_find == std::get<0>(array[guess]))
		{

			//posto moze biti vise istih minimizatora (samo na razlicitim pozicijama) moramo se vratiti natrag do prvoga
			result = guess;

			while (result > 0 && std::get<0>(array[result - 1]) == number_to_find)
				result--;

			return result;
		}
		else if (std::get<0>(array[guess]) < number_to_find) {
			min = guess + 1;
		}
		else {
			max = guess - 1;
		}
	}

	return result;
}

void minimizer_matches(
	std::unordered_map <unsigned int, unsigned int> &reference_genome_minimizers_positions_in_vector,
        std::vector < std::vector <std::tuple <unsigned int, unsigned int, bool>>> &reference_genome_minimizers,
	std::vector <std::tuple <unsigned int, unsigned int, bool>> &sequence_minimizers,
	std::vector <std::tuple <unsigned int, unsigned int>> &matches_same_strand,
	std::vector <std::tuple <unsigned int, unsigned int>> &matches_different_strand
) {

	for (auto tuple : sequence_minimizers) {
		if (std::get<2>(tuple) == true) {
			if (reference_genome_minimizers_positions_in_vector.find(std::get<0>(tuple)) != reference_genome_minimizers_positions_in_vector.end())
				for (auto tupleR : reference_genome_minimizers[reference_genome_minimizers_positions_in_vector[std::get<0>(tuple)] - 1]) {
					if (std::get<2>(tupleR) == true)
						matches_same_strand.emplace_back(std::make_tuple(std::get<1>(tupleR), std::get<1>(tuple)));
					else
						matches_different_strand.emplace_back(std::make_tuple(std::get<1>(tupleR), std::get<1>(tuple)));
				}
		}
		else {
			if (reference_genome_minimizers_positions_in_vector.find(std::get<0>(tuple)) != reference_genome_minimizers_positions_in_vector.end())
				for (auto tupleR : reference_genome_minimizers[reference_genome_minimizers_positions_in_vector[std::get<0>(tuple)] - 1]) {
					if (std::get<2>(tupleR) == true)
                                                matches_different_strand.emplace_back(std::make_tuple(std::get<1>(tupleR), std::get<1>(tuple)));
                                        else
                                                matches_same_strand.emplace_back(std::make_tuple(std::get<1>(tupleR), std::get<1>(tuple)));
				}
		}

	}

}

bool compare_by_diagonal(std::tuple<unsigned int, unsigned int> a, std::tuple<unsigned int, unsigned int> b) {
	return ((int)std::get<1>(a) - (int)std::get<0>(a)) < ((int)std::get<1>(b) - (int)std::get<0>(b));
}

bool compare_q_t(std::tuple<unsigned int, unsigned int> a, std::tuple<unsigned int, unsigned int> b) {
	if (std::get<0>(a) != std::get<0>(b))
		return std::get<0>(a) < std::get<0>(b);
	else
		return std::get<1>(a) < std::get<1>(b);
}

std::vector <std::vector <std::tuple <unsigned int, unsigned int>>> prepare_match_groups(std::vector <std::tuple <unsigned int, unsigned int>> &matches/*, unsigned int &no_matches*/) {
	std::sort(matches.begin(), matches.end(), compare_by_diagonal);
	std::vector <std::tuple <unsigned int, unsigned int>> match_group;
	std::vector <std::vector <std::tuple <unsigned int, unsigned int>>> match_groups;

	unsigned int size = matches.size();
	int last_diagonal;
	int next_diagonal;
	unsigned int i = 0;
	unsigned int match_group_size = 0;

	if (size == 0) {
		//std::cerr << "\nWarning! No minimizer matches.\n";
		//no_matches++;
		return match_groups;
	}

	do {

		match_group.emplace_back(matches[i]);
		match_group_size = 1;
		last_diagonal = std::get<1>(matches[i]) - std::get<0>(matches[i]);
		i++;

		if (i < size) {
			for (next_diagonal = std::get<1>(matches[i]) - std::get<0>(matches[i]); last_diagonal + 500 > next_diagonal; next_diagonal = std::get<1>(matches[i]) - std::get<0>(matches[i])) {
				match_group.emplace_back(matches[i]);
				match_group_size++;
				last_diagonal = next_diagonal;
				if (++i == size)
					break;
			}
		}

		if (match_group_size > 4) {
			std::sort(match_group.begin(), match_group.end(), compare_q_t);
			match_groups.emplace_back(match_group);
		}

		match_group.clear();

	} while (i < size);

	return match_groups;

}

//funkcija koja vraca indeks elementa koji je najmanji element veci ili jednak key (naravno funkcija funkcionira jedino ako je vektor
//sortiran i provjerili smo da je key nije veci od najveceg elementa vektora (tj. ne vrijedi: key > najveci element))
unsigned int CeilIndex(std::vector <std::tuple <unsigned int, unsigned int, unsigned int>> &v, int l, int r, unsigned int key) {
	while (r - l > 1) {
		int m = l + (r - l) / 2;

		if (std::get<2>(v[m]) >= key)
			r = m;
		else
			l = m;
	}

	return r;
}

//funkcija koja trazi LIS za zadane matcheve minimizatora (funkcija pazi da u LIS-u ne budu matchevi s istom vrijednoscu query -> pomocno polje)
std::pair<int, int> find_LIS(std::vector <std::tuple <unsigned int, unsigned int>> &minimizer_matches) {
	std::vector <std::tuple <unsigned int, unsigned int, unsigned int>> tail(minimizer_matches.size()); //polje za biljezit LIS-a
	//polja su redom:
	//	1. pozicija prvog matcha LIS-a u minimizer_matches
	//	2. pozicija zadnjeg matcha LIS-a u minimizer_matches
	//	3. zadnji element LIS-a

	std::vector <std::tuple <unsigned int, unsigned int, unsigned int>> auxiliary_tail(minimizer_matches.size()); //pomocno polje, ovo polje sluzi kako bi eliminirali
														//pojavljivanje matcheva koji imaju istu vrijednost query

	std::vector <unsigned int> changed_indexes; //vektor koji biljezi gdje su se dogodile promjene u tail tako da ih mozemo prepisati u auxiliary_tail
						//umjesto uvijek da kad predemo sve matcheve koji imaju isti query_position prepisemo cijeli tail u
						//auxiliary_tail mozemo samo prepisati promijenjene elemente
	if (minimizer_matches.size() == 0)
		return std::make_pair(-1, -1); //ako nema matcheva ne treba radit alignment

	unsigned int length = 1; //broji duljinu LIS-a
	unsigned int auxiliary_tail_length = 1; //velicina pomocnog polja

	//minimizer_matches su sortirani po query_position(prva vrijednost) pa po target_position(druga vrijednost)
	//kao prvi element LIS-a postavljamo match koji ima najmanju vrijednost target_position, ali
	//gledamo samo one elemente s prvom (najmanjom) vrijednost query_position

	unsigned int query_value_min = std::get<0>(minimizer_matches[0]);
	unsigned int first_target_value_position = 0;

	for (unsigned int i = 1; query_value_min == std::get<0>(minimizer_matches[i]); i++)
		if (std::get<1>(minimizer_matches[first_target_value_position]) < std::get<1>(minimizer_matches[i]))
			first_target_value_position = i;

	tail[0] = std::make_tuple(first_target_value_position, first_target_value_position, std::get<1>(minimizer_matches[first_target_value_position]));
	auxiliary_tail[0] = tail[0];

	//	for (int i = 0; i < length; i++)
	//		std::cout << std::get<0>(tail[i]) << std::get<1>(tail[i]) << std::get<2>(tail[i]) << std::endl;

	//std::cout << std::endl;

	//sad prolazimo kroz ostale matcheve i racunamo LIS
	//da bi eliminari pojavu elemenata s istom vrijednosti query koristimo pomocno polje

	for (unsigned int i = 1; i < minimizer_matches.size(); i++) {

		unsigned int target_position = std::get<1>(minimizer_matches[i]);

		//manji od najmanjeg -> samo zamjenimo najmanji element
		if (target_position <= std::get<2>(tail[0])) {
			tail[0] = std::make_tuple(i, i, target_position);
			changed_indexes.emplace_back(0);
		}
		//veci od najveceg: dvije opcije
		//1 -> ako najveci element nema isti query_position kao i novi element, onda dodamo novi element u vektor
		//2 -> ako najveci element ima isti query_position kao i novi element, onda provjerimo postoji li element na istoj poziciji 
		//	-> kao i najveci element u pomocnom vektoru:
		//	-> ako da, onda stvorimo novi element u tail ako je on veci od elemnta u pomocnom vektoru
		//  -> ako ne, onda ne dodajemo nista
		else if (target_position > std::get<2>(tail[length - 1])) {
			if (std::get<0>(minimizer_matches[std::get<1>(tail[length - 1])]) != std::get<0>(minimizer_matches[i])) {
				tail[length] = std::make_tuple(std::get<0>(tail[length - 1]), i, target_position);
				changed_indexes.emplace_back(length);
				length++;
			}
			else if (length == auxiliary_tail_length && std::get<2>(auxiliary_tail[length - 1]) < target_position) {
				tail[length] = std::make_tuple(std::get<0>(auxiliary_tail[length - 1]), i, target_position);
				changed_indexes.emplace_back(length);
				length++;
			}
		}

		//ako nije nista od prethodnog, onda trazimo indeks najmanjeg veceg ili jednakog broja od novog elementa -> dobijemo indeks
		//promatramo element na indeks - 1
		//dvije opcije:
		//1 -> ako taj element nema isti query_position kao i novi element, onda zamijenimo element na indeksu: indeks sa novim elementom
		//2 -> ako taj element ima isti query_position kao i novi element, onda pogledamo u pomocnom vektoru element na poziciji indeks - 1
		//dvije opcije:
		//1 -> ako je taj element manji od novog elementa, ako zamijenimo u vektoru tail element na indeksu: indeks s novim elementom
		//2 -> ako ne, ne radimo nista

		else {
			unsigned int index = CeilIndex(tail, -1, length - 1, target_position);

			if (std::get<0>(minimizer_matches[std::get<1>(tail[index - 1])]) != std::get<0>(minimizer_matches[i])) {
				tail[index] = std::make_tuple(std::get<0>(tail[index - 1]), i, target_position);
				changed_indexes.emplace_back(index);
			}

			else if (std::get<2>(auxiliary_tail[index - 1]) < target_position) {
				tail[index] = std::make_tuple(std::get<0>(auxiliary_tail[index - 1]), i, target_position);
				changed_indexes.emplace_back(index);
			}

		}

		//nakon sto prodemo sve elemente s istim query_position moramo prekopirati tail u auxiliary_tail
		if (i + 1 < minimizer_matches.size() && std::get<0>(minimizer_matches[i]) != std::get<0>(minimizer_matches[i + 1])) {
			for (auto index : changed_indexes)
				auxiliary_tail[index] = tail[index];

			changed_indexes.clear();

			auxiliary_tail_length = length;
		}


		/*		for (int i = 0; i < length; i++)
					std::cout << std::get<0>(tail[i]) << std::get<1>(tail[i]) << std::get<2>(tail[i]) << std::endl;

				std::cout << std::endl;

				for (int i = 0; i < auxiliary_tail_length; i++)
					std::cout << std::get<0>(auxiliary_tail[i]) << std::get<1>(auxiliary_tail[i]) << std::get<2>(auxiliary_tail[i]) << std::endl;

				std::cout << std::endl;
		*/
	}

	if (length >= 4)
		return std::make_pair(std::get<0>(tail[length - 1]), std::get<1>(tail[length - 1]));
	else
		return std::make_pair(-1, -1);
}


void minimizer_occurrences (std::vector<std::unique_ptr<SequenceFormat>> &sequences, unsigned int k, unsigned int window_size, float f) {
	std::unordered_map <unsigned int, unsigned int> minimizer_occurrences;
	std::vector <std::tuple<unsigned int, unsigned int, bool>> current_minimizers;

	for (auto &ptr : sequences) {
		current_minimizers = white::minimizers(ptr -> sequence.c_str(), ptr -> sequence.size(), k, window_size);

		for (auto minimizer_tuple : current_minimizers) {
			minimizer_occurrences[std::get<0>(minimizer_tuple)]++;
		}
	}

	/*std::ofstream fout;
  	fout.open ("minimizer_occurrences.csv", std::ios::out);

	if (!fout.is_open()){
		std::cout << "Unable to open file.\n";
		exit(1);
	}

  	fout << "Minimizer,Number of occurrences\r\n"; //posto mi je linux na windowsima da mogu citat i u notepadu

	for (std::unordered_map<unsigned int, unsigned int>::iterator it = minimizer_occurrences.begin(); it != minimizer_occurrences.end(); it++)
	{
		if(it->second != 1) {
  			fout << it->first << ",";
			fout << it->second << "\r\n";
		}
	}

  	fout.close();*/

	unsigned int number_of_minimizers_to_disregard = (unsigned int) std::round(minimizer_occurrences.size() * f);
	std::vector<unsigned int> frequencies;
	unsigned int number_of_minimizers_with_frequency_of_1 = 0;

	for (std::unordered_map<unsigned int, unsigned int>::iterator it = minimizer_occurrences.begin(); it != minimizer_occurrences.end(); it++) {
		frequencies.emplace_back (it -> second);

		if (it -> second == 1)
			number_of_minimizers_with_frequency_of_1++;
	}

	std::sort (frequencies.begin(), frequencies.end());

	printf("\n--- Minimizer statistics ---\n\n"

		"-> Minimizers found: %lu									\n"
		"-> Number of minimizers with frequency of 1: %u						\n"
		"-> After disregarding top %f minimizers the frequency of the most frequent minimizer is:	\n"
		"	Frequency: %u										\n\n",

		minimizer_occurrences.size(), number_of_minimizers_with_frequency_of_1, f, frequencies[frequencies.size() - 1 - number_of_minimizers_to_disregard]);
}

//pomocna funkcija koja obavlja racunanje LIS-a nad match group, radi alignment i vraca PAF koji treba ispisati u obliku stringa
std::string/* void */ LIS_alignment_PAF(std::vector <std::tuple <unsigned int, unsigned int>> &match_group,
						std::unique_ptr<SequenceFormat> &seq,
						std::unique_ptr<SequenceFormat> &ref,
						char strand,
						int match,
						int mismatch,
						int gap,
						int k,
						bool print_cigar
) {

	std::pair <long long int, long long int> LIS_begin_end; //ovo su pozicije prvog i zadnjeg minimizer match-a LIS-a u match_group-u
	std::string query;
	unsigned int q_begin;
	unsigned int q_end;
	unsigned int q_length;

	std::string target;
	unsigned int t_begin;
	unsigned int t_end;
	unsigned int t_length;

	LIS_begin_end = find_LIS(match_group);

	if (std::get<0>(LIS_begin_end) == -1 && std::get<1>(LIS_begin_end) == -1)
		return "Empty";

	q_begin = std::get<1>(match_group[std::get<0>(LIS_begin_end)]);
	q_end = std::get<1>(match_group[std::get<1>(LIS_begin_end)]) + (k - 1); //za jedan i drugi end treba dodat jos "k-1" jer se
										//moramo pozicionirat na kraj zadnjeg minimizatora
	t_begin = std::get<0>(match_group[std::get<0>(LIS_begin_end)]);
	t_end = std::get<0>(match_group[std::get<1>(LIS_begin_end)]) + (k - 1);

	std::ostringstream os;
	std::string paf_sekvence;

	if (print_cigar) {

		std::string cigar;
        	unsigned int target_begin = 0;

		q_length = q_end - q_begin + 1;
		t_length = t_end - t_begin + 1;

		query = (seq->sequence).substr(q_begin, q_length);
		target = (ref->sequence).substr(t_begin, t_length);

		int value = white::pairwise_alignment(query.c_str(), q_length, target.c_str(), t_length,
			white::AlignmentType::global, match, mismatch, gap, cigar, target_begin);

		//stvaranje PAF

		//izracun broja matcheva u alignmentu

		int cigar_size = cigar.size();
		int number_of_matches = 0;

		for (int i = 0; i < cigar_size; i++)
			if (cigar[i] == '=')
				number_of_matches++;

		os << seq->name.c_str() << "\t" << seq->sequence.size() << "\t" << q_begin << "\t" << q_end << "\t" << strand << "\t" << ref -> name.c_str() << "\t" << ref -> sequence.size();
		os << "\t" << t_begin << "\t" << t_end << "\t" << number_of_matches << "\t" << cigar_size << "\t" << "255";

		paf_sekvence = os.str();

		//ispis PAF
/*		printf("\n%s\t%lu\t%d\t%d\t%c\t%s\t%lu\t%d\t%d\t%d\t%d\t%d",
			seq->name.c_str(),
			seq->sequence.size(),
			q_begin,
			q_end,
			strand,
			ref -> name.c_str(),
			ref -> sequence.size(),
			t_begin,
			t_end,
			number_of_matches,
			cigar_size,
			255);

		printf("\tcg:Z:<%s>\n", cigar.c_str());
*/	}

	else {

		os << seq->name.c_str() << "\t" << seq->sequence.size() << "\t" << q_begin << "\t" << q_end << "\t" << strand << "\t" << ref -> name.c_str() << "\t" << ref -> sequence.size();
                os << "\t" << t_begin << "\t" << t_end << "\t" << "0" << "\t" << "0" << "\t" << "255";

                paf_sekvence = os.str();

		//ispis PAF
/*               printf("\n%s\t%lu\t%d\t%d\t%c\t%s\t%lu\t%d\t%d\t%d\t%d\t%d",
                        seq->name.c_str(),
                        seq->sequence.size(),
                        q_begin,
                        q_end,
                        strand,
                        ref -> name.c_str(),
                        ref -> sequence.size(),
                        t_begin,
                        t_end,
                        0,
                        0,
                        255);

                printf("\n");*/
	}

	return paf_sekvence;
}

//pomocna funkcija kako bi mogli napraviti efikasniji paralelizam

void map_sequence_to_reference(
	std::unordered_map <unsigned int, unsigned int> &reference_genome_minimizers_positions_in_vector,
        std::vector < std::vector <std::tuple <unsigned int, unsigned int, bool>>> &reference_genome_minimizers,
	std::vector <std::unique_ptr <SequenceFormat>> &sequences,
	std::vector <std::unique_ptr <SequenceFormat>> &reference_genome,
	int match,
	int mismatch,
	int gap,
	int k,
	int window_size,
	bool print_cigar,
	unsigned int start,
	unsigned int end,
	std::vector <std::string> &ispis_za_pojedinu_dretvu
) {

	//unsigned int no_matches = 0;
	//unsigned int no_match_groups_original = 0;
	//unsigned int no_match_groups_reverse_complement = 0;

	for (unsigned int i = start; i < end; i++) {
		std::vector <std::tuple <unsigned int, unsigned int, bool>> sequence_minimizer_index;

        	std::vector <std::tuple <unsigned int, unsigned int>> matches_same_strand;
		std::vector <std::tuple <unsigned int, unsigned int>> matches_different_strand;

        	std::vector <std::vector <std::tuple <unsigned int, unsigned int>>> match_groups_same_strand;
	        std::vector <std::vector <std::tuple <unsigned int, unsigned int>>> match_groups_different_strand;

		//za sekvencu racunamo minimizatore
		sequence_minimizer_index = white::minimizers(sequences[i]->sequence.c_str(), sequences[i]->sequence.size(), k, window_size);

		//radimo matcheve minimizatora (i za minimizatore koji dolaze iz originalnih sekvenci i reverznih komplementa)
		//te ih spremamo u matches_original i matches_reverse_complement

		minimizer_matches(
			reference_genome_minimizers_positions_in_vector,
			reference_genome_minimizers,
			sequence_minimizer_index,
			matches_same_strand,
			matches_different_strand
		);

		//prepare match_groups sortira matcheve najprije po dijagonali, zatim radi grupe nad kojima cemo traziti LIS (te grupe ce biti
		//sortirane po poziciji untar target (reference) i poziciji unutar query (sort za matcheve s jednakim target position))
		match_groups_same_strand = prepare_match_groups(matches_same_strand);
		match_groups_different_strand = prepare_match_groups(matches_different_strand);

		//za svaku match grupu racunamo LIS, zatim radimo alignment dobivenih regija i printamo PAF
		if (!match_groups_same_strand.empty())
			for (auto match_group : match_groups_same_strand)
				ispis_za_pojedinu_dretvu.emplace_back (LIS_alignment_PAF (match_group, sequences[i], reference_genome[0], '+', match, mismatch, gap, k, print_cigar));
		//else
		//	no_match_groups_original++;

		if(!match_groups_different_strand.empty())
			for (auto match_group : match_groups_different_strand)
				ispis_za_pojedinu_dretvu.emplace_back (LIS_alignment_PAF (match_group, sequences[i], reference_genome[0], '-', match, mismatch, gap, k, print_cigar));
		//else
		//	no_match_groups_reverse_complement++;
	}
		//std::cout << no_matches << std::endl;
		//std::cout << no_match_groups_original << " " << no_match_groups_reverse_complement << std::endl;
}

//funkcija koja razdvaja minimizatore reference u one s originalnog stranda i reverznog stranda

void map_reference_genome_minimizers (
	std::vector <std::tuple <unsigned int, unsigned int, bool>> &reference_minimizer_index,
	std::unordered_map <unsigned int, unsigned int> &reference_genome_minimizers_positions_in_vector,
        std::vector < std::vector <std::tuple <unsigned int, unsigned int, bool>>> &reference_genome_minimizers
){

	unsigned int i = 0;

        std::vector <std::tuple <unsigned int, unsigned int, bool>> temp;

        for (auto tuple : reference_minimizer_index) {
		if (reference_genome_minimizers_positions_in_vector[std::get<0>(tuple)] == 0) {
                	i++;
                        temp.emplace_back(tuple);
                        reference_genome_minimizers_positions_in_vector[std::get<0>(tuple)] = i;
                        reference_genome_minimizers.emplace_back(temp);
                        temp.clear();
                }
                else
                        reference_genome_minimizers[reference_genome_minimizers_positions_in_vector[std::get<0>(tuple)] - 1].emplace_back(tuple);
	}
}


int main (int argc, char* argv[])
{
	int option = 0;
	std::vector<std::unique_ptr<SequenceFormat>> sequences;
	std::vector<std::unique_ptr<SequenceFormat>> reference_genome;

	int match = MATCH;
        int mismatch = MISMATCH;
        int gap = GAP;
	int k = K_MER_LENGTH;
	int window_size = WINDOW_SIZE;
	float f = PERCENTAGE;

	bool print_cigar = false;
	unsigned int number_of_threads = NUMBER_OF_THREADS;

        while ((option = getopt_long(argc, argv, "hvm:s:g:k:w:f:ct:", long_options, NULL)) != -1) {
                switch (option) {
                        case 'h':
				help();
                                exit(0);

                        case 'v':
				version();
                                exit(0);

			case 'm':
				match = atoi (optarg);
				break;

			case 's':
                                mismatch = atoi (optarg);
				break;

			case 'g':
                                gap = atoi (optarg);
				break;

			case 'k':
				k = atoi (optarg);
				break;

			case 'w':
				window_size = atoi (optarg);
				break;

			case 'f':
				f = std::stof (optarg, NULL);
				break;

			case 'c':
				print_cigar = true;
				break;

			case 't':
				number_of_threads = atoi(optarg);
				break;

                        default:
                                fprintf (stderr, "\n--Error:\n\tUnsupported option. Usage: %s [-h] [--help] [-v] [--version] [-m ARG] [-s ARG] [-g ARG] [-k ARG] [-w ARG] [-f ARG] "
						"sequence_file reference_genome_file match mismatch gap\n\n", argv[0]);
                                exit(0);

                }
	}

	if (argc - optind == 2)
	{
		if (check_format (argv[optind], fasta_formats))
			sequences = parse_file<bioparser::FastaParser> (std::string (argv[optind]));

		else if (check_format (argv[optind], fastq_formats))
                        sequences = parse_file<bioparser::FastqParser> (std::string (argv[optind]));

		else
			fprintf (stderr, "\nError:													\n"
					 "	Unsuported format. File containing sequences (1st file) needs to have one of the following extensions:	\n"
                                     	 "	-> .fasta             -> .fastq										\n"
                                     	 "	-> .fa                -> .fq										\n"
                                     	 "	-> .fasta.gz          -> .fastq.gz									\n"
                                     	 "	-> .fa.gz             -> .fq.gz										\n\n");

		if (check_format (argv[++optind], fasta_formats))
                        reference_genome = parse_file<bioparser::FastaParser> (std::string (argv[optind]));
		else

			fprintf (stderr, "\n --Error:														\n"
					 "	Unsuported format. File containing reference genome (2nd file) needs to have one of the following extensions:	\n"
                                     	 "	-> .fasta													\n"
                                     	 "	-> .fa														\n"
                                     	 "	-> .fasta.gz													\n"
                                     	 "	-> .fa.gz													\n\n");
	}

	 else
                        fprintf (stderr, "\n--Error:\n"
					 "	Usage: white_mapper [OPTIONS] sequence_file reference_genome_file match mismatch gap	\n\n"

					 "	Number of nonoption arguments must be 2:						\n"
                                         "		1st file contains sequences in FASTA or FASTQ format				\n"
                                         "		2nd file contains a reference genome in FASTA fromat				\n\n");

	//to use these next few lines change constructor of SequenceFormat to public and minimizer_occurrences (sekvenceTest, k, window_size)
	/*std::vector<std::unique_ptr<SequenceFormat>> sekvenceTest;
	std::unique_ptr<SequenceFormat> p1 (new SequenceFormat ("S1", 2, "TCAGGAAGAAGCAGA", 15));
	std::unique_ptr<SequenceFormat> p2 (new SequenceFormat ("S2", 2, "GTCATGCACGTTCAC", 15));
	std::unique_ptr<SequenceFormat> p3 (new SequenceFormat ("S3", 2, "TCAGGAAGAAGCAGA", 15));
	sekvenceTest.push_back(std::move(p1));
	sekvenceTest.push_back(std::move(p2));
	sekvenceTest.push_back(std::move(p3));*/

	//finding minimizers and making a csv file of their occurrences
/*	minimizer_occurrences (sequences, k, window_size, f);

	//alignment of 2 random selected sequences
	srand (time (NULL));

	int i1 = rand() % sequences.size();
	int i2 = rand() % sequences.size();

	//printf("Sequences used:\n%d_%d", i1, i2);

	int values[3];
*/
//	std::string cigar;
//	unsigned int target_begin = 0;
/*	//printf("\nFirst sequence (Length: %lu):\n%s\n\nSecond sequence (Length: %lu):\n%s\n\n", sequences[i1] -> sequence.size(), sequences[i1] -> sequence.c_str(),
	//					sequences[i2] -> sequence.size(), sequences[i2] -> sequence.c_str());

	int value = white::pairwise_alignment(	sequences[i1] -> sequence.c_str(), sequences[i1] -> sequence.size(),
  	                                   	sequences[i2] -> sequence.c_str(), sequences[i2] -> sequence.size(),
 	                                  	white::AlignmentType::global, match, mismatch, gap, cigar, target_begin);

	values[0] = value;

	std::cout << "\nGlobal alignment:\n";
	//std::cout << value << std::endl;
	//std::cout << target_begin << std::endl;
	printf ("%s\n", cigar.c_str());

	value = white::pairwise_alignment(  sequences[i1] -> sequence.c_str(), sequences[i1] -> sequence.size(),
                                            sequences[i2] -> sequence.c_str(), sequences[i2] -> sequence.size(),
                                            white::AlignmentType::semi_global, match, mismatch, gap, cigar, target_begin);

	values[1] = value;

        std::cout << "\nSemi-global alignment:\n";
        //std::cout << value << std::endl;
        //std::cout << target_begin << std::endl;
	printf ("%s\n", cigar.c_str());

	value = white::pairwise_alignment(  sequences[i1] -> sequence.c_str(), sequences[i1] -> sequence.size(),
                                            sequences[i2] -> sequence.c_str(), sequences[i2] -> sequence.size(),
                                            white::AlignmentType::local, match, mismatch, gap, cigar, target_begin);

	values[2] = value;

        std::cout << "\nLocal alignment:\n";
        //std::cout << value << std::endl;
        //std::cout << target_begin << std::endl;
	printf ("%s\n", cigar.c_str());
*/

/*	std::string q = {"TCCG"};
	std::string t = {"ACTCCGAT"};
	std::string cigar;
	unsigned int target_begin = 0;
	int value = white::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), white::AlignmentType::semi_global, match, mismatch, gap, cigar, target_begin);
	std::cout << value << std::endl;
	std::cout << target_begin << std::endl;
	std::cout << cigar << std::endl;
	std::cout <<std::endl;
*/
/*	std::string q = {"TTCCGCCAA"};
	std::string t = {"AACCCCTT"};
	std::string cigar;
	unsigned int target_begin = 0;
	int value = white::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), white::AlignmentType::local, match, mismatch, gap, cigar, target_begin);
	std::cout << value << std::endl;
	std::cout << target_begin << std::endl;
	std::cout << cigar << std::endl;
	std::cout <<std::endl;
*/

/*	std::string q = {"TGCATAT"};
	std::string t = {"ATCCGAT"};
	std::string cigar;
	unsigned int target_begin = 0;
	int value = white::pairwise_alignment(q.c_str(), q.size(), t.c_str(), t.size(), white::AlignmentType::global, match, mismatch, gap, cigar, target_begin);
	std::cout << value << std::endl;
	std::cout << target_begin << std::endl;
	std::cout << cigar << std::endl;
	std::cout <<std::endl;
*/


//-----------------------------------------------ZADNJI_ZADATAK-----------------------------------------------------------------------------------------------------------------------

	std::vector <std::tuple <unsigned int, unsigned int, bool>> reference_minimizer_index = reference_genome_minimizers(reference_genome[0], k, window_size, f);
	std::vector <std::tuple <unsigned int, unsigned int, bool>> sequence_minimizer_index;
//---------
	std::unordered_map <unsigned int, unsigned int> reference_genome_minimizers_positions_in_vector;
        std::vector < std::vector <std::tuple <unsigned int, unsigned int, bool>>> reference_genome_minimizers;

	map_reference_genome_minimizers (
		reference_minimizer_index,
		reference_genome_minimizers_positions_in_vector,
		reference_genome_minimizers
	);

//---------

//	for (auto tuple : reference_minimizer_index)
//		std::cout << std::get<0>(tuple) << " " << std::get<1>(tuple) << " " << (std::get<2>(tuple) ? "true" : "false") << std::endl;

	//std::string reference_name = reference_genome[0] -> name;
	//int reference_size = reference_genome[0] -> sequence.size();

	//stvorimo thread_pool kako bi paralelizirali izvodenje koda
	std::shared_ptr < thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool();

	std::vector<std::future<void>> thread_futures;

	unsigned int size = sequences.size();

	unsigned int sequences_per_thread = size / number_of_threads;
	unsigned int residue_sequences = size % number_of_threads;

	unsigned int start = 0;
	unsigned int end = 0;

	//vector za pohranjivati PAF-ove koje na kraju izrsavanja svih dretvi treba isprintati
	std::vector <std::vector <std::string>> ispis;

	for (unsigned int i = 0; i < number_of_threads; i++) {
		std::vector <std::string> ispis_za_pojedinu_dretvu;
		ispis.emplace_back(ispis_za_pojedinu_dretvu);
	}

	for (unsigned int i = 0; i < number_of_threads; i++) {
		if (residue_sequences > 0) {
			end = start + sequences_per_thread + 1;
			residue_sequences--;
		}
		else
			end = start + sequences_per_thread;

		thread_futures.emplace_back(thread_pool->submit_task(
			map_sequence_to_reference,
			std::ref(reference_genome_minimizers_positions_in_vector),
			std::ref(reference_genome_minimizers),
			std::ref(sequences),
			std::ref(reference_genome),
			match,
			mismatch,
			gap,
			k,
			window_size,
			print_cigar,
			start,
			end,
			std::ref(ispis[i])
		));

		start = end;
	}

	for (auto &it : thread_futures)
		it.wait();

	for (auto ispis_za_pojedinu_dretvu : ispis)
		for (auto paf_sekvence : ispis_za_pojedinu_dretvu){
			if (paf_sekvence.compare("Empty") != 0)
				printf("%s\n", paf_sekvence.c_str());
		}

	return 0;
}
