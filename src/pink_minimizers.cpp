#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <unordered_set>
#include <queue>

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
    
    
    void get_interior_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length,std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers_set){
        
        unsigned int l = window_length + k - 1;
        unsigned int last_beginning = sequence_length - l;
        
        unsigned int first_min_value;
        unsigned int first_min_value_reversed;
        unsigned int temp_value;
        unsigned int temp_value_reversed;
        
        std::deque<std::tuple<unsigned int, unsigned int, bool>> red;
        
        first_min_value = get_first_value(sequence, 0, k);
        first_min_value_reversed = get_first_value_reversed(sequence, 0, k);
        
        if(first_min_value < first_min_value_reversed) {
            red.emplace_back(first_min_value, 0, 0);
        } else {
            red.emplace_back(first_min_value_reversed, 0, 1);
        }
        temp_value = first_min_value;
        temp_value_reversed = first_min_value_reversed;
        
        for (unsigned int i = 1; i < window_length-1; i++) {
            temp_value = get_value(sequence, i, k, temp_value);
            temp_value_reversed = get_value_reversed(sequence, i, k, temp_value_reversed);
            if(temp_value < temp_value_reversed) {
                red.emplace_back(temp_value, i, 0);
            } else {
                red.emplace_back(temp_value_reversed, i, 1);
            }
        }
        
        
        
        for (unsigned beginning_position = window_length-1; beginning_position <= last_beginning; beginning_position++){
            
            temp_value = get_value(sequence, beginning_position, k, temp_value);
            temp_value_reversed = get_value_reversed(sequence, beginning_position, k, temp_value_reversed);
            int min_orientation_tmp = temp_value < temp_value_reversed ? 0 : 1;
            int min_value_tmp = temp_value < temp_value_reversed ? temp_value : temp_value_reversed;
            
            std::tuple<unsigned int, unsigned int, bool> tmp_elem = std::make_tuple(min_value_tmp, beginning_position, min_orientation_tmp);
            
            while(!red.empty()) {
                std::tuple<unsigned int, unsigned int, bool> current_elemn = red.back();
                if(std::get<0>(current_elemn) >= std::get<0>(tmp_elem)) {
                    red.pop_back();
                } else {
                    break;
                }
            }
            
            red.push_back(tmp_elem);
            
            
            std::tuple<unsigned int, unsigned int, bool> tmp_rez = red.front();
            minimizers_set.insert(tmp_rez);
            
            if(std::get<1>(tmp_rez) + (window_length-1) <= beginning_position) {
                red.pop_front();
            }
        }
    }
    
    void get_beginning_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers_set){

        unsigned int temp_value;
        unsigned int temp_value_reversed;
        
        unsigned int first_min_value = get_first_value(sequence, 0, k);
        unsigned int first_min_value_reversed = get_first_value_reversed(sequence, 0, k);
        
        std::deque<std::tuple<unsigned int, unsigned int, bool>> red;
        
        
        if(first_min_value < first_min_value_reversed) {
            red.emplace_back(first_min_value, 0, 0);
        } else {
            red.emplace_back(first_min_value_reversed, 0, 1);
        }
        temp_value = first_min_value;
        temp_value_reversed = first_min_value_reversed;
        
        for (unsigned int i = 1; i < window_length-1; i++) {
            temp_value = get_value(sequence, i, k, temp_value);
            temp_value_reversed = get_value_reversed(sequence, i, k, temp_value_reversed);
            if(temp_value < temp_value_reversed) {
                red.emplace_back(temp_value, i, 0);
            } else {
                red.emplace_back(temp_value_reversed, i, 1);
            }
        }
        
        std::tuple<unsigned int, unsigned int, bool> tmp_rez = red.front();
        minimizers_set.insert(tmp_rez);
        
    }
    
    void get_end_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers_set){
      
        unsigned int last_beginning = sequence_length - k;
        unsigned int last_window = sequence_length - window_length;
        
        unsigned int first_min_value = get_first_value(sequence, last_beginning, k);
        unsigned int first_min_value_reversed = get_first_value_reversed(sequence, last_beginning, k);
        unsigned int temp_value;
        unsigned int temp_value_reversed;
        
        std::deque<std::tuple<unsigned int, unsigned int, bool>> red;
        
        
        if(first_min_value < first_min_value_reversed) {
            red.emplace_back(first_min_value, last_beginning, 0);
        } else {
            red.emplace_back(first_min_value_reversed, last_beginning, 1);
        }
        temp_value = first_min_value;
        temp_value_reversed = first_min_value_reversed;
        
        for (unsigned int i = last_beginning; i >= last_window; i--) {
            temp_value = get_value(sequence, i, k, temp_value);
            temp_value_reversed = get_value_reversed(sequence, i, k, temp_value_reversed);
            
            if(temp_value < temp_value_reversed) {
                red.emplace_back(temp_value, i, 0);
            } else {
                red.emplace_back(temp_value_reversed, i, 1);
            }
        }
        
        std::tuple<unsigned int, unsigned int, bool> tmp_rez = red.front();
        minimizers_set.insert(tmp_rez);
        
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


