#pragma once

namespace brown {

std::vector<std::pair<unsigned int, unsigned int>> minimizers(
    const char* sequence,
    unsigned int sequence_length,
    unsigned int k,
    unsigned int window_length);

}