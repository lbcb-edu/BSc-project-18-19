#ifndef BLUE_MINIMIZERS_HPP_INCLUDED
#define BLUE_MINIMIZERS_HPP_INCLUDED

#include <vector>
#include <tuple>

namespace blue
{
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(
                const char* sequence, unsigned int sequence_length,
                unsigned int k, unsigned int window_length);
}

#endif // BLUE_MINIMIZERS_HPP_INCLUDED
