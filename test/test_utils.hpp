
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>

std::string to_hex( std::string const& str ) ;

template<typename ... Args>
std::string format( std::string const& format, Args... args )
{
	size_t size = std::snprintf( nullptr, 0, format.c_str(), args... ) + 1; // Extra space for '\0'
	std::unique_ptr<char[]> buf( new char[ size ] ) ; 
	std::snprintf( buf.get(), size, format.c_str(), args... );
	return std::string( buf.get(), buf.get() + size - 1 );
}

