#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <unordered_set>

#include "pink_minimizers.hpp"

namespace std {
    template <> struct hash<tuple<unsigned int, unsigned int, bool >> {
        inline size_t operator()(const tuple<unsigned int, unsigned int, bool > &value) const {
            return hash<int>()(get<0>(value) ^ get<1>(value) ^ get<2>(value));
        }
    };
}
namespace pink {
    unsigned int get_letter_value(char letter){
        switch (letter) {
            case 'C':
                return 0;
            case 'A':
                return 1;
            case 'T':
                return 2;
            case 'G':
                return 3;
            default:
                return 4;
        }
    }
    unsigned int get_complement_letter_value(char letter){
        switch (letter) {
            case 'C':
                return 3;
            case 'A':
                return 2;
            case 'T':
                return 1;
            case 'G':
                return 0;
            default:
                return 4;
        }
    }
    unsigned int get_first_value(const char* sequence, unsigned int position, unsigned int k){
        unsigned int value = 0b0;
        
        for (unsigned i = position; i < position + k; i++){
            value = value | get_letter_value(sequence[i]);
            if (i < position + k - 1)
                value = value << 2;
        }
        return value;
    }
    
    unsigned int get_value(const char* sequence, unsigned int position, unsigned int k, unsigned int temp_value){
        unsigned int value = 0b0;
        unsigned int mask = (1 << (2*k)) - 1;
        
        value = ((temp_value << 2) & mask) | get_letter_value(sequence[position + k - 1]);
        
        return value;
    }
    
    unsigned int get_first_value_reversed(const char* sequence, unsigned int position, unsigned int k){
        unsigned int value = 0b0;
        for (unsigned i = position + k; i > position; i--){
            value = value | get_complement_letter_value(sequence[i - 1]);
            if (i > position + 1)
                value = value << 2;
        }
        return value;
    }
        
    
    unsigned int get_value_reversed(const char* sequence, unsigned int position, unsigned int k, unsigned int temp_value_reversed){
        unsigned int value = 0b0;
            value = ((get_complement_letter_value(sequence[position+k-1]) << 2*k) | temp_value_reversed) >> 2;
        return value;
    }
    
    void calculate_min_value(unsigned int* min_value, unsigned int* min_position, bool* min_direction, const char* sequence, unsigned int minimizer_beginning, unsigned int k,unsigned int* temp_value, unsigned int* temp_value_reversed){
        
        *temp_value = get_value(sequence, minimizer_beginning, k, *temp_value);
        *temp_value_reversed = get_value_reversed(sequence, minimizer_beginning, k, *temp_value_reversed);

        if (*temp_value < *min_value){
            *min_value = *temp_value;
            *min_position = minimizer_beginning;
            *min_direction = 0;
        }
        if (*temp_value_reversed < *min_value){
            *min_value = *temp_value_reversed;
            *min_position = minimizer_beginning;
            *min_direction = 1;
        }
    }
    
    void get_interior_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length,std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers_set){
        
        unsigned int l = window_length + k - 1;
        unsigned int last_beginning = sequence_length - l;
        
        unsigned int min_value;
        unsigned int min_position;
        bool min_direction;
        unsigned int minimizer_beginning;
        
        unsigned int first_min_value = get_first_value(sequence, 0, k);
        unsigned int first_min_value_reversed = get_first_value_reversed(sequence, 0, k);
        
        for (unsigned beginning_position = 0; beginning_position <= last_beginning; beginning_position++){
            
            
            
            unsigned int temp_value = 0;
            unsigned int temp_value_reversed = 0;
            
            for (unsigned i = 0; i < window_length; i++){
                minimizer_beginning = beginning_position + i;
                if (i == 0){
                    if (first_min_value < first_min_value_reversed){
                        min_value = first_min_value;
                        min_position = beginning_position;
                        min_direction = 0;
                    } else {
                        min_value = first_min_value_reversed;
                        min_position = beginning_position;
                        min_direction = 1;
                    }
                    
                    temp_value = first_min_value;
                    temp_value_reversed = first_min_value_reversed;
                    
                } else {
                    calculate_min_value(&min_value, &min_position, &min_direction, sequence, minimizer_beginning, k, &temp_value, &temp_value_reversed);
                }
                if (i == 1){
                    first_min_value = temp_value;
                    first_min_value_reversed = temp_value_reversed;
                }
                
            }
            minimizers_set.emplace(min_value, min_position, min_direction);
        }
    }
    
    void get_beginning_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers_set){
        
        unsigned int min_value;
        unsigned int min_position;
        bool min_direction;
        unsigned int minimizer_beginning;
        
        unsigned int first_min_value = get_first_value(sequence, 0, k);
        unsigned int first_min_value_reversed = get_first_value_reversed(sequence, 0, k);
        
        for (unsigned beginning_position = 1; beginning_position < window_length - 1; beginning_position++){
            min_value = get_value(sequence, 0, k, 0);
            min_position = 0;
            min_direction = 0;
            unsigned int temp_value = 0;
            unsigned int temp_value_reversed = 0;
            
            for (unsigned i = 0; i < beginning_position; i++){
                if (i == 0){
                    
                    if (first_min_value < first_min_value_reversed){
                        min_value = first_min_value;
                        min_position = beginning_position - 1;
                        min_direction = 0;
                    } else {
                        min_value = first_min_value_reversed;
                        min_position = beginning_position - 1;
                        min_direction = 1;
                    }
                    
                    temp_value = first_min_value;
                    temp_value_reversed = first_min_value_reversed;
                } else {
                    minimizer_beginning = beginning_position + i;
                    calculate_min_value(&min_value, &min_position, &min_direction, sequence, minimizer_beginning, k, &temp_value, &temp_value_reversed);
                }
                if (i == 1){
                    first_min_value = temp_value;
                    first_min_value_reversed = temp_value_reversed;
                }
            }
            minimizers_set.emplace(min_value, min_position, min_direction);
        }
    }
    
    void get_end_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers_set){
        
        unsigned int min_value;
        unsigned int min_position;
        bool min_direction;
        unsigned int minimizer_beginning;
        unsigned int last_beginning = sequence_length - k;
        unsigned int last_window = sequence_length - window_length;
        
        unsigned int first_min_value = get_first_value(sequence, last_beginning - 1, k);
        unsigned int first_min_value_reversed = get_first_value_reversed(sequence, last_beginning - 1, k);
        
        for (unsigned int beginning_position = last_beginning - 1; beginning_position > last_window; beginning_position--){
            min_value = get_value(sequence, last_beginning, k, last_beginning);
            min_position = last_beginning;
            min_direction = 0;
            unsigned int temp_value = 0;
            unsigned int temp_value_reversed = 0;
            
            
            for (unsigned i = beginning_position; i < last_beginning; i++){
                if (i == last_beginning - 1){
                    
                    if (first_min_value < first_min_value_reversed){
                        min_value = first_min_value;
                        min_position = beginning_position;
                        min_direction = 0;
                    } else {
                        min_value = first_min_value_reversed;
                        min_position = beginning_position;
                        min_direction = 1;
                    }
                    temp_value = first_min_value;
                    temp_value_reversed = first_min_value_reversed;
                } else {
                    minimizer_beginning = beginning_position + i;
                    calculate_min_value(&min_value, &min_position, &min_direction, sequence, minimizer_beginning, k, &temp_value, &temp_value_reversed);
                }
                if (i == last_beginning){
                    first_min_value = temp_value;
                    first_min_value_reversed = temp_value_reversed;
                }
            }
            minimizers_set.emplace(min_value, min_position, min_direction);
        }
    }
    
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length){
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector;
        std::unordered_set<std::tuple<unsigned int, unsigned int, bool>> minimizers_set;
        
        get_interior_minimizers(sequence, sequence_length, k, window_length, minimizers_set);
        get_beginning_minimizers(sequence, sequence_length, k, window_length, minimizers_set);
        get_end_minimizers(sequence, sequence_length, k, window_length, minimizers_set);
        
        for (auto minimizer: minimizers_set)
            minimizers_vector.emplace_back(minimizer);
        
        return minimizers_vector;
        
    }
}
