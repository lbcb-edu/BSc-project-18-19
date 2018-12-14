#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <utility>
#include <tuple>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <map>

namespace blue
{

    long long int getMask(int length);
    unsigned long long int calculateFirstWindow(int length, const char* sequence);
    unsigned long long int calculateReverseWindow(int length, const char* sequence, unsigned int sequence_length);
    std::set<std::tuple<unsigned int, unsigned int, bool>> findMinimizers(unsigned long long int kMers, long long int mask,
                                                                          unsigned int k, unsigned int window_length,
                                                                          int position, bool isOriginal);

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(
                    const char* sequence, unsigned int sequence_length,
                    unsigned int k, unsigned int window_length) {

            std::set<std::tuple<unsigned int, unsigned int, bool>> uniqueMinimizers;

            int length = window_length + k - 1;
            int endCond = sequence_length - length;
            long long int windowMask = getMask(length);
            long long int mask = getMask(k);
            const char* temp = sequence + sequence_length - 1;
            unsigned long long int endWindow = calculateReverseWindow(length, temp, sequence_length);
            unsigned long long int window = calculateFirstWindow(length, sequence);
            sequence = sequence + length;
            temp = temp - length;

            for(int i = 0; i <= endCond; i++) {
                //find minimizers from original
                std::set<std::tuple<unsigned int, unsigned int, bool>> tempOriginal = findMinimizers(window, mask, k, window_length, i, true);

                //find minimizers from reverse complement
                std::set<std::tuple<unsigned int, unsigned int, bool>> tempRCompl = findMinimizers(endWindow, mask, k, window_length, i, false);

                if(std::get<0>(*tempOriginal.begin()) > std::get<0>(*tempRCompl.begin())) {
                    uniqueMinimizers.insert(tempRCompl.begin(), tempRCompl.end());

                } else if(std::get<0>(*tempOriginal.begin()) < std::get<0>(*tempRCompl.begin())) {
                    uniqueMinimizers.insert(tempOriginal.begin(), tempOriginal.end());

                } else {
                    uniqueMinimizers.insert(tempOriginal.begin(), tempOriginal.end());
                    uniqueMinimizers.insert(tempRCompl.begin(), tempRCompl.end());
                }

                unsigned int code;
                switch(*sequence) {
                    case 'C':
                        code = 0;
                        break;
                    case 'A':
                        code = 1;
                        break;
                    case 'T':
                        code = 2;
                        break;
                    case 'G':
                        code = 3;
                        break;
                    default:
                        break;
                }

                unsigned int rCode = 0;
                switch(*temp) {
                    case 'C':
                        rCode = 3;
                        break;
                    case 'A':
                        rCode = 2;
                        break;
                    case 'T':
                        rCode = 1;
                        break;
                    case 'G':
                        rCode = 0;
                        break;
                    default:
                        break;
                }

                endWindow = ((endWindow << 2) | rCode) & windowMask;
                --temp;
                window = ((window << 2) | code) & windowMask;
                ++sequence;
            }

            std::vector<std::tuple<unsigned int, unsigned int, bool>> output(uniqueMinimizers.begin(), uniqueMinimizers.end());
            return output;
    }

    //returns set od minimizers for given string
    std::set<std::tuple<unsigned int, unsigned int, bool>> findMinimizers(unsigned long long int kMers, long long int mask, unsigned int k, unsigned int window_length, int position, bool isOriginal) {

        std::set<std::tuple<unsigned int, unsigned int, bool>> minimizerSet;
        //set first kmer as minimum
        unsigned int minimum = kMers & mask;
        minimizerSet.emplace(minimum, position+window_length-1, isOriginal);

        for(int i = 1; i <= window_length; i++) {
            unsigned int bitResult = kMers & mask;
            kMers = kMers >> 2;
            if(bitResult > minimum) continue;
            if(bitResult < minimum) {
                minimum = bitResult;
                minimizerSet.clear();
            }

            minimizerSet.emplace(bitResult, position+window_length-i, isOriginal);
        }
        return minimizerSet;
    }

    unsigned long long int calculateFirstWindow(int length, const char* sequence) {
        unsigned long long int bitResult = 0;

        for(int i = 0; i < length; i++) {
        char c = *sequence;
        unsigned int digit;
            switch(c) {
                case 'C':
                    digit = 0;
                    break;
                case 'A':
                    digit = 1;
                    break;
                case 'T':
                    digit = 2;
                    break;
                case 'G':
                    digit = 3;
                    break;
                default:
                    std::cout << "Sequence is not valid!" << std::endl;
            }
            c=*++sequence;
            bitResult = bitResult << 2;
            bitResult = bitResult | digit;
        }
        return bitResult;
    }

    unsigned long long int calculateReverseWindow(int length, const char* temp, unsigned int sequence_length) {
        unsigned long long int bitResult = 0;
        char c = *temp;

        for(int i = 0; i < length; i++) {
            unsigned int digit;
            switch(c) {
                case 'C':
                    digit = 3;
                    break;
                case 'A':
                    digit = 2;
                    break;
                case 'T':
                    digit = 1;
                    break;
                case 'G':
                    digit = 0;
                    break;
                default:
                    std::cout << "Sequence is not valid!" << std::endl;
            }
            c=*--temp;
            bitResult = bitResult << 2;
            bitResult = bitResult | digit;
        }
        return bitResult;
    }

    long long int getMask(int length) {
        long long int mask = 0;
        for(int i = 0; i < length; i++) {
            mask = mask << 2;
            mask = mask | 3;
        }
        return mask;
    }
}
