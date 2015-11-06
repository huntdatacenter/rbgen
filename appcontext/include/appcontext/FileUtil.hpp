
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef APPCONTEXT_FILEUTIL_HPP
#define APPCONTEXT_FILEUTIL_HPP

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <memory>
#include <iterator>
#include <iostream>

namespace appcontext {
	typedef std::auto_ptr< std::istream > INPUT_FILE_PTR ;
	typedef std::auto_ptr< std::ostream > OUTPUT_FILE_PTR ;

	std::string determine_file_compression( std::string const& filename ) ;
	std::string determine_file_mode( std::string const& filename ) ;

	// Return a stream representing a given output file, optionally with compression.
	OUTPUT_FILE_PTR
	open_file_for_output( std::string const& filename, std::string const& compression, std::string const& mode = "text" ) ;

	// Return a stream representing a given output file, guessing compression option from the filename.
	OUTPUT_FILE_PTR
	open_file_for_output( std::string const& filename, std::string const& mode = "text" ) ;
}

#endif

