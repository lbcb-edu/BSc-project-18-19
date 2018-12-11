#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <tuple>

#include "pink_minimizers.hpp"

namespace pink {
    
    unsigned int get_value(const char* sequence, unsigned int position, unsigned int k, std::map<char, unsigned int> odd_position){
        unsigned int value = 0;
        for (unsigned i = position; i < position + k; i++){
	    //std::cout << sequence[i] << std::endl;
            if (i % 2 == 0)
                value = 10*value + (3 - odd_position[sequence[i]]);
            else
                value = 10*value + odd_position[sequence[i]];
        }
	//std::cout << "value " << value << std::endl;
        return value;
    }
    
    unsigned int get_value_reversed(const char* sequence, unsigned int sequence_length, unsigned int position, unsigned int k, std::map<char, unsigned int> odd_position){
        unsigned int value = 0;

		/*std::cout << sequence << std::endl;
		std::cout << sequence_length << std::endl;
		std::cout << position << std::endl;
		std::cout << k << std::endl;*/

        for (unsigned i = position + k - 1; i > position - 1; i--){

		//std::cout << sequence[i] << std::endl;
            if ((sequence_length + 1 - i) % 2 == 0)
                value = 10*value + odd_position[sequence[i]]; 
            else
                value = 10*value + (3 - odd_position[sequence[i]]);
        }
	//std::cout << "value r " << value << std::endl;
        return value;
    }
    
    void calculate_min_value(unsigned int* min_value, unsigned int* min_position, bool* min_direction, const char* sequence,unsigned int sequence_length, unsigned int minimizer_beginning, unsigned int k, std::map<char, unsigned int> odd_position){

	//std::cout << minimizer_beginning << std::endl;
        
        unsigned int temp_value = get_value(sequence, minimizer_beginning, k, odd_position);
        unsigned int temp_value_reversed = get_value_reversed(sequence, sequence_length, minimizer_beginning, k, odd_position);

	
        
        if(temp_value < *min_value){
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
    
    std::vector<std::tuple<unsigned int, unsigned int, bool>> get_interior_minimizers(const char* sequence, unsigned int sequence_length,
                                                                                       unsigned int k,
                                                                                      unsigned int window_length, std::map<char, unsigned int> odd_position){
        std::tuple<unsigned int, unsigned int, bool> minimizers_tuple;
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector;
        
        unsigned int l = window_length + k - 1;
        unsigned int last_beginning = sequence_length - l;
        
        unsigned int min_value;
        unsigned int min_position;
        bool min_direction;
        unsigned int minimizer_beginning;
        
        for (unsigned beginning_position = 0; beginning_position < last_beginning; beginning_position++){
            

            min_value = get_value(sequence, beginning_position, k,odd_position);
            min_position = beginning_position;
            min_direction = 0;
            
            for (unsigned i = 0; i <= window_length; i++){
                
                minimizer_beginning = beginning_position + i;
              //  std::cout << minimizer_beginning << std::endl;
                calculate_min_value(&min_value, &min_position, &min_direction, sequence, sequence_length, minimizer_beginning, k, odd_position);
                
            }
            minimizers_tuple = std::make_tuple(min_value, min_position, min_direction);
	    std::cout << "tuple " << min_value << " " << min_position << " " << min_direction << std::endl;
            minimizers_vector.push_back(minimizers_tuple);
        }
        
        return minimizers_vector;
        
    }
    
    std::vector<std::tuple<unsigned int, unsigned int, bool>> get_beginning(const char* sequence, unsigned int sequence_length,
                                                                            unsigned int k,
                                                                            unsigned int window_length, std::map<char, unsigned int> odd_position){
        
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector;
        std::tuple<unsigned int, unsigned int, bool> minimizers_tuple;
        
        unsigned int last_beginning = sequence_length - k;
        unsigned int min_value;
        unsigned int min_position;
        bool min_direction;
        unsigned int minimizer_beginning;
        
        for (unsigned beginning_position = 1; beginning_position < last_beginning; beginning_position++){
            
            min_value = get_value(sequence, 0, k ,odd_position);
            min_position = 0;
            min_direction = 0;
            
            for (unsigned i = 0; i <= beginning_position; i++){
                
                minimizer_beginning = beginning_position + i;
                
                calculate_min_value(&min_value, &min_position, &min_direction, sequence,sequence_length, minimizer_beginning, k, odd_position);
                
            }
            minimizers_tuple = std::make_tuple(min_value, min_position, min_direction);
	    std::cout << "tuple " << min_value << " " << min_position << " " << min_direction << std::endl;
            minimizers_vector.push_back(minimizers_tuple);
        }
        
        return minimizers_vector;
    }
    
    
    std::vector<std::tuple<unsigned int, unsigned int, bool>> get_end(const char* sequence, unsigned int sequence_length,
                                                                            unsigned int k,
                                                                      unsigned int window_length, std::map<char, unsigned int> odd_position){
        
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector;
        std::tuple<unsigned int, unsigned int, bool> minimizers_tuple;
        
        unsigned int min_value;
        unsigned int min_position;
        bool min_direction;
        unsigned int minimizer_beginning;
        unsigned int last_beginning = sequence_length - k;
        
        for (unsigned beginning_position = last_beginning - 1; beginning_position >= 0; beginning_position--){
            
            min_value = get_value(sequence, last_beginning, k, odd_position);
            min_position = last_beginning;
            min_direction = 0;
            
            for (unsigned i = beginning_position; i <= last_beginning; i++){
                
                minimizer_beginning = beginning_position + i;
                
                calculate_min_value(&min_value, &min_position, &min_direction, sequence, sequence_length, minimizer_beginning, k, odd_position);
                
            }
            minimizers_tuple = std::make_tuple(min_value, min_position, min_direction);
	    std::cout << "tuple " << min_value << " " << min_position << " " << min_direction << std::endl;
            minimizers_vector.push_back(minimizers_tuple);
        }
        
        return minimizers_vector;
    }
    
    std::vector<std::tuple<unsigned int, unsigned int, bool>> get_end_minimizers(const char* sequence, unsigned int sequence_length,
                                                                                      unsigned int k,
                                                                                 unsigned int window_length, std::map<char, unsigned int> odd_position){
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector;
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector_beginning;
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector_end;
        
        minimizers_vector_beginning = get_beginning(sequence, sequence_length, k, window_length, odd_position);
        minimizers_vector_end = get_end(sequence, sequence_length, k, window_length, odd_position);
        
        minimizers_vector = minimizers_vector_beginning;
        minimizers_vector.insert(minimizers_vector.end(), minimizers_vector_end.begin(), minimizers_vector_end.end());
        
        return minimizers_vector;
        
    }
    
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length,
                                                                         unsigned int k,
                                                                         unsigned int window_length){
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector;
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector_interior;
        std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector_end;
        std::map<char, unsigned int> odd_position;
        
        odd_position['C'] = 0;
        odd_position['A'] = 1;
        odd_position['T'] = 2;
        odd_position['U'] = 2;
        odd_position['G'] = 3;
        
        
        minimizers_vector_interior = get_interior_minimizers(sequence, sequence_length, k, window_length, odd_position);
        minimizers_vector_end = get_end_minimizers(sequence, sequence_length, k, window_length, odd_position);
        
        minimizers_vector = minimizers_vector_interior;
        minimizers_vector.insert(minimizers_vector.end(), minimizers_vector_end.begin(), minimizers_vector_end.end());
        
        return minimizers_vector;
        
    }
    
}
/*int main() {
	
    const char* sequence = "CGAC";
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers_vector;
    
    minimizers_vector = pink::minimizers(sequence, 4, 3, 1); //?
    
    for (auto const& minimizer: minimizers_vector){
        std::cout << std::get<0>(minimizer) << " " << std::get<1>(minimizer) << " " << std::get<2>(minimizer) << std::endl;
    }
    return 0;
}*/
