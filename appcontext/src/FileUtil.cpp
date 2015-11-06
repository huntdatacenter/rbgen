
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <cassert>
#include <stdexcept>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
namespace bio = boost::iostreams ;

#include <boost/filesystem.hpp>
namespace BFS = boost::filesystem ;
#include "appcontext/FileUtil.hpp"

namespace appcontext {
	// Return a stream representing a given output file, optionally with gzip compression.
	OUTPUT_FILE_PTR
	open_file_for_output( std::string const& filename, std::string const& compression, std::string const& mode ) {
		std::ios_base::openmode open_mode = std::ios_base::out ;
		if( mode == "binary" ) {
			open_mode |= std::ios_base::binary ;
		}
		std::auto_ptr< bio::filtering_ostream > stream_ptr( new bio::filtering_ostream ) ;
		if( compression == "gzip" ) {
			stream_ptr->push( bio::gzip_compressor() ) ;
		} else if( compression == "bzip2" ) {
			assert(0) ; // not implemented yet
		} else if( compression == "no compression" ){
			// nothing to do
		} else {
			assert(0) ;
		}
		bio::file_sink sink( filename, open_mode ) ;
		if( !sink.is_open() ) {
			throw std::runtime_error( "appcontext::open_file_for_output(): could not open file \"" + filename + "\"." ) ;
		}
		stream_ptr->push( sink ) ;

		return OUTPUT_FILE_PTR( stream_ptr.release() ) ;
	}

	std::string determine_file_compression( std::string const& filename ) {
		if( filename.find( ".gz" ) == ( filename.size()-3 )) {	
			return "gzip" ;
		}
		else if( filename.find( ".bz2" ) == ( filename.size()-4 )) {
			return "bzip2" ;
		}
		else {
			return "no compression" ;
		}
	}

	// Return a stream representing a given input file, attempting to auto-detect the compression to
	// use from the filename.
	OUTPUT_FILE_PTR
	open_file_for_output( std::string const& filename, std::string const& mode ) {
		return open_file_for_output( filename, determine_file_compression( filename ), mode ) ;
	}
}
