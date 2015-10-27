
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "test_utils.hpp"

std::string to_hex( std::string const& str ) {
	std::ostringstream o ;
	for( std::size_t i = 0; i < str.size(); ++i ) {
		if( i % 4 == 0 )
			o << "|" ;
		o << std::hex << std::setw(2) << std::setfill('0') << static_cast<int> ( static_cast<unsigned char>( str[i] ) ) ;
	}
	return o.str() ;
}

