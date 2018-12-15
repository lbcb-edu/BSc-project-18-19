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
            return hash<int>()(get<0>(value) ^ get<0>(value) ^ get<0>(value));
        }
    };
}

namespace pink {
    
unsigned int get_value(const char* sequence, unsigned int position, unsigned int k, std::unordered_map<char, unsigned int> position_value){
    unsigned int value = 0b0;

    for (unsigned i = position; i < position + k; i++){
        value = value | position_value[sequence[i]];
        if (i < position + k - 1)
            value = value << 2;
    }
    return value;
}

unsigned int get_value_reversed(const char* sequence, unsigned int sequence_length, unsigned int position, unsigned int k, std::unordered_map<char, unsigned int> position_value){
    unsigned int value = 0b0;
    
    for (unsigned i = position + k; i > position; i--){
        value = value | (3 - position_value[sequence[i - 1]]);
        if (i > position + 1)
            value = value << 2;
    }
    return value;
}

void calculate_min_value(unsigned int* min_value, unsigned int* min_position, bool* min_direction, const char* sequence,unsigned int sequence_length, unsigned int minimizer_beginning, unsigned int k, std::unordered_map<char, unsigned int> position_value){
    
    unsigned int temp_value = get_value(sequence, minimizer_beginning, k, position_value);
    unsigned int temp_value_reversed = get_value_reversed(sequence, sequence_length, minimizer_beginning, k, position_value);
    
    if (temp_value < *min_value){
        *min_value = temp_value;
        *min_position = minimizer_beginning;
        *min_direction = 0;
    }

    if (temp_value_reversed < *min_value){
        *min_value = temp_value_reversed;
        *min_position = minimizer_beginning;
        *min_direction = 1;
    }  
}
    
void get_interior_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length, std::unordered_map<char, unsigned int> position_value,std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers_set){
    
    unsigned int l = window_length + k - 1;
    unsigned int last_beginning = sequence_length - l;
    
    unsigned int min_value;
    unsigned int min_position;
    bool min_direction;
    unsigned int minimizer_beginning;
    
    for (unsigned beginning_position = 0; beginning_position <= last_beginning; beginning_position++){
        
        
        min_value = get_value(sequence, beginning_position, k, position_value);
        min_position = beginning_position;
        min_direction = 0;
        
        for (unsigned i = 0; i < window_length; i++){
            
            minimizer_beginning = beginning_position + i;
            calculate_min_value(&min_value, &min_position, &min_direction, sequence, sequence_length, minimizer_beginning, k, position_value);
            
        }
        minimizers_set.emplace(min_value, min_position, min_direction);
    }
    
}
    
void get_beginning_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length, std::unordered_map<char, unsigned int> position_value, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers_set){
    
    unsigned int min_value;
    unsigned int min_position;
    bool min_direction;
    unsigned int minimizer_beginning;
    
    for (unsigned beginning_position = 1; beginning_position < window_length - 1; beginning_position++){
        
        min_value = get_value(sequence, 0, k , position_value);
        min_position = 0;
        min_direction = 0;
        
        for (unsigned i = 0; i < beginning_position; i++){
            
            minimizer_beginning = beginning_position + i;
            
            calculate_min_value(&min_value, &min_position, &min_direction, sequence,sequence_length, minimizer_beginning, k, position_value);
        }
        minimizers_set.emplace(min_value, min_position, min_direction);
        
    }
}
    
void get_end_minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length, std::unordered_map<char, unsigned int> position_value, std::unordered_set<std::tuple<unsigned int, unsigned int, bool>>& minimizers_set){
    
    unsigned int min_value;
    unsigned int min_position;
    bool min_direction;
    unsigned int minimizer_beginning;
    unsigned int last_beginning = sequence_length - k;
    unsigned int last_window = sequence_length - window_length;
    
    for (unsigned int beginning_position = last_beginning - 1; beginning_position > last_window; beginning_position--){
        
        min_value = get_value(sequence, last_beginning, k, position_value);
        min_position = last_beginning;
        min_direction = 0;
        
        for (unsigned i = beginning_position; i < last_beginning; i++){
            
            minimizer_beginning = beginning_position + i;
            
            calculate_min_value(&min_value, &min_position, &min_direction, sequence, sequence_length, minimizer_beginning, k, position_value);
            
        }
        minimizers_set.emplace(min_value, min_position, min_direction);
    }
}
    
std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length,
                                                                     unsigned int k,
                                                                     unsigned int window_length){
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector;
    std::unordered_set<std::tuple<unsigned int, unsigned int, bool>> minimizers_set;
    std::unordered_map<char, unsigned int> position_value;
    
    position_value['C'] = 0;
    position_value['A'] = 1;
    position_value['T'] = 2;
    position_value['U'] = 2;
    position_value['G'] = 3;
    
    
    get_interior_minimizers(sequence, sequence_length, k, window_length, position_value, minimizers_set);
    get_beginning_minimizers(sequence, sequence_length, k, window_length, position_value, minimizers_set);
    get_end_minimizers(sequence, sequence_length, k, window_length, position_value, minimizers_set);
    
    for (auto minimizer: minimizers_set)
        minimizers_vector.emplace_back(minimizer);
    
    return minimizers_vector;
    
}
}
//int main() {
//    const char* sequence = "TGACGTACATGGACA";
//    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector;
//
//    minimizers_vector = pink::minimizers(sequence, 15, 3, 3);
//
//    for (auto const& minimizer: minimizers_vector){
//        std::cout << std::get<0>(minimizer) << " " << std::get<1>(minimizer) << " " << std::get<2>(minimizer) << std::endl;
//    }
//    return 0;
//}

