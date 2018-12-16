#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <utility>
#include <tuple>
#include <algorithm>
#include <unordered_set>

namespace std {
template <> struct hash<std::tuple<unsigned int, unsigned int, bool >> {
    inline size_t operator()(const std::tuple<unsigned int, unsigned int, bool > &v) const {
        std::hash<int> int_hasher;
        return int_hasher(std::get<0>(v)) ^ int_hasher(std::get<1>(v)) ^ int_hasher(std::get<2>(v));
    }
};
}

namespace blue
{
    long long int getMask(int length);
    unsigned long long int calculateFirstWindow(int length, const char* sequence);
    unsigned long long int calculateReverseWindow(int length, const char* sequence, unsigned int sequence_length);
    std::vector<std::tuple<unsigned int, unsigned int, bool>> findMinimizers(unsigned long long int kMers, long long int mask,
                                                                          unsigned int k, unsigned int window_length,
                                                                          int position, bool isOriginal, std::tuple<unsigned int, unsigned int, bool> current);

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(
                    const char* sequence, unsigned int sequence_length,
                    unsigned int k, unsigned int window_length) {

            std::unordered_set<std::tuple<unsigned int, unsigned int, bool>> uniqueMinimizers;

            int length = window_length + k - 1;
            int endCond = sequence_length - length;
            long long int windowMask = getMask(length);
            long long int mask = getMask(k);
            const char* temp = sequence + sequence_length - 1;
            unsigned long long int endWindow = calculateReverseWindow(length, temp, sequence_length);
            unsigned long long int window = calculateFirstWindow(length, sequence);
            sequence = sequence + length;
            temp = temp - length;
            std::tuple<unsigned int, unsigned int, bool> currentO = std::make_tuple(-1, -1, 0);
            std::tuple<unsigned int, unsigned int, bool> currentR = std::make_tuple(-1, -1, 0);
            uniqueMinimizers.emplace((window & (mask << 2*(window_length - 1))) >> 2*(window_length - 1), 0, 1);

            for(int i = 0; i <= endCond; i++) {
                //find minimizers from original
                std::vector<std::tuple<unsigned int, unsigned int, bool>> tempOriginal = findMinimizers(window, mask, k, window_length, i, true, currentO);

                //find minimizers from reverse complement
                std::vector<std::tuple<unsigned int, unsigned int, bool>> tempRCompl = findMinimizers(endWindow, mask, k, window_length, i, false, currentR);

                if(std::get<0>(tempOriginal[0]) > std::get<0>(tempRCompl[0])) {
                    uniqueMinimizers.insert(tempRCompl.begin(), tempRCompl.end());
                    currentR = tempRCompl[0];

                } else if(std::get<0>(tempOriginal[0]) < std::get<0>(tempRCompl[0])) {
                    uniqueMinimizers.insert(tempOriginal.begin(), tempOriginal.end());
                    currentO = tempOriginal[0];

                } else {
                    uniqueMinimizers.insert(tempOriginal.begin(), tempOriginal.end());
                    uniqueMinimizers.insert(tempRCompl.begin(), tempRCompl.end());
                    currentO = tempOriginal[0];
                    currentR = tempRCompl[0];
                }

                if (i == endCond) {
                    break;
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

            uniqueMinimizers.emplace(window & mask, sequence_length - k, 1);

            std::vector<std::tuple<unsigned int, unsigned int, bool>> output (uniqueMinimizers.begin(), uniqueMinimizers.end());
            return output;
    }

    //returns set od minimizers for given string
    std::vector<std::tuple<unsigned int, unsigned int, bool>> findMinimizers(unsigned long long int kMers, long long int mask, unsigned int k, unsigned int window_length, int position,
     bool isOriginal, std::tuple<unsigned int, unsigned int, bool> current) {

        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizerSet;
        //set first kmer as minimum
        unsigned int minimum = kMers & mask;

        if (std::get<1>(current) >= position) {
            if (std::get<0>(current) < minimum) {
                minimizerSet.emplace_back(current);
                return minimizerSet;
            }
            if (std::get<0>(current) == minimum) {
                minimizerSet.emplace_back(minimum, position+window_length-1, isOriginal);
                return minimizerSet;
            }
        }

        minimizerSet.emplace_back(minimum, position+window_length-1, isOriginal);

        for(int i = 1; i <= window_length; i++) {
            unsigned int bitResult = kMers & mask;
            kMers = kMers >> 2;
            if(bitResult > minimum) continue;
            if(bitResult < minimum) {
                minimum = bitResult;
                minimizerSet.clear();
            }

            minimizerSet.emplace_back(bitResult, position+window_length-i, isOriginal);
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
