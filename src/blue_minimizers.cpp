#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <utility>
#include <tuple>
#include <algorithm>

namespace blue
{
    unsigned int calculateBits(std::string kmer, unsigned int k);
    std::set<std::tuple<unsigned int, unsigned int, bool>> findMinimizers(std::string kMers, unsigned int k, int position, bool isOriginal);
    std::string makeNumberString(int flag, const char* sequence);
    std::string reverseComplement(std::string original);

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(
                    const char* sequence, unsigned int sequence_length,
                    unsigned int k, unsigned int window_length) {

            std::string numberString = makeNumberString(k%2, sequence);
            std::set<std::tuple<unsigned int, unsigned int, bool>> uniqueMinimizers;

            int length = window_length + k - 1;
            int endCond = sequence_length - length;
            for(int i = 0; i <= endCond; i++) {
                //find minimizers from original
                std::string original = numberString.substr(i, length);
                std::set<std::tuple<unsigned int, unsigned int, bool>> tempOriginal = findMinimizers(original, k, i, true);

                //find minimizers from reverse complement
                std::string reverseCompl = reverseComplement(original);
                std::set<std::tuple<unsigned int, unsigned int, bool>> tempRCompl = findMinimizers(reverseCompl, k, i, false);

                if(std::get<0>(*tempOriginal.begin()) > std::get<0>(*tempRCompl.begin())) {
                    uniqueMinimizers.insert(tempRCompl.begin(), tempRCompl.end());

                } else if(std::get<0>(*tempOriginal.begin()) < std::get<0>(*tempRCompl.begin())) {
                    uniqueMinimizers.insert(tempOriginal.begin(), tempOriginal.end());

                } else {
                    uniqueMinimizers.insert(tempOriginal.begin(), tempOriginal.end());
                    uniqueMinimizers.insert(tempRCompl.begin(), tempRCompl.end());
                }
            }

            //end-minimizers
            unsigned int startEnd = calculateBits(numberString.substr(0, k), k);
            unsigned int endEnd = calculateBits(numberString.substr(sequence_length-k, k), k);

            uniqueMinimizers.insert(std::make_tuple(startEnd, 0, true));
            uniqueMinimizers.insert(std::make_tuple(endEnd, sequence_length-k, true));

            std::vector<std::tuple<unsigned int, unsigned int, bool>> output(uniqueMinimizers.begin(), uniqueMinimizers.end());
            return output;
        }

    //returns coded string representation of DNA/RNA sequence
    std::string makeNumberString(int flag, const char* sequence) {
        std::string numberString = "";
        for(char c = *sequence; c != '\0'; c=*++sequence) {
            switch(c) {
                case 'C':
                    numberString += flag == 1 ? "0" : "3";
                    break;
                case 'A':
                    numberString += flag == 1 ? "1" : "2";
                    break;
                case 'T':
                    numberString += flag == 1 ? "2" : "1";
                    break;
                case 'G':
                    numberString += flag == 1 ? "3" : "0";
                    break;
                default:
                    std::cout << "Sequence is not valid!" << std::endl;
            }

        }
        return numberString;
    }

    //returns set od minimizers for given string
    std::set<std::tuple<unsigned int, unsigned int, bool>> findMinimizers(std::string kMers, unsigned int k, int position, bool isOriginal) {

        std::set<std::tuple<unsigned int, unsigned int, bool>> minimizerSet;
        //set first kmer as minimum
        unsigned int minimum = calculateBits(kMers.substr(0, k), k);
        minimizerSet.insert(std::make_tuple(minimum, position, isOriginal));

        for(int i = 1, length = kMers.length() - k; i <= length; i++) {
            std::string kmer = kMers.substr(i, k);
            unsigned int bitResult = calculateBits(kmer, k);

            if(bitResult > minimum) continue;
            if(bitResult < minimum) {
                minimum = bitResult;
                minimizerSet.clear();
            }

            minimizerSet.insert(std::make_tuple(bitResult, position+i, isOriginal));
        }

        if(position == 0)
        return minimizerSet;
    }

    //returns integer representations of kmer
    unsigned int calculateBits(std::string kmer, unsigned int k) {
        unsigned int bitResult = 0;

        for(int j = 0; j < k; j++) {
            unsigned int digit = (unsigned int)kmer.at(j) - '0';
            bitResult = bitResult << 2;
            bitResult = bitResult | digit;
        }

        return bitResult;
    }

    //returns reversed complement of DNA/RNA string sequence
    std::string reverseComplement(std::string original) {
        std::string numberString = "";
        for(int i = 0, length = original.length(); i < length; i++) {
            switch(original.at(i)) {
                case '0':
                    numberString += "3";
                    break;
                case '1':
                    numberString += "2";
                    break;
                case '2':
                    numberString += "1";
                    break;
                case '3':
                    numberString += "0";
                    break;
                default:
                    std::cout << "Sequence is not valid!" << std::endl;
            }
        }
        reverse(numberString.begin(), numberString.end());
        return numberString;
    }
}
