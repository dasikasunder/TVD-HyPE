/*
 * utilities.cc
 *      Author: sunder
 */

#include "headers.hh"

namespace Utilities {
std::string int_to_string (unsigned int value, const unsigned int digits) {
    std::string lc_string = std::to_string(value);

    if (lc_string.size() < digits) {
        // We have to add the padding zeroes in front of the number
        const unsigned int padding_position = (lc_string[0] == '-')
                                                ?
                                                1
                                                :
                                                0;

        const std::string padding(digits - lc_string.size(), '0');
        lc_string.insert(padding_position, padding);
    }

    return lc_string;
}
}


