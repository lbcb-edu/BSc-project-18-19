#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <algorithm>
#include <list>
#include <utility>

namespace orange {

    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(const char* sequence, unsigned int sequence_length, unsigned int k, unsigned int window_length);

}