#pragma once

namespace white {

	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers(	const char* sequence,
										unsigned int sequence_length,
                                                                     		unsigned int k,
                                                                     		unsigned int window_length);
}
