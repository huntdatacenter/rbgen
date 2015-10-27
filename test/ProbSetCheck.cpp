
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include "catch.hpp"
#include "genfile/types.hpp"
#include "test_utils.hpp"
#include "ProbSetCheck.hpp"

/*
* ProbSetCheck
* This class checks sequence and values of probabilities set using
* the bgen API
*/
ProbSetCheck::ProbSetCheck(
	std::size_t n,
	GetExpectedProbs get_expected_probs
):
	m_number_of_samples( n ),
	m_get_expected_probs( get_expected_probs ),
	m_sample_i( std::numeric_limits< std::size_t >::max() ),
	m_number_of_entries( std::numeric_limits< std::size_t >::max() ),
	m_entry_i( std::numeric_limits< std::size_t >::max() ),
	m_state( eNone )
{}

ProbSetCheck::~ProbSetCheck() throw() {
	REQUIRE( m_sample_i + 1 == m_number_of_samples ) ;
	REQUIRE( m_entry_i == m_number_of_entries ) ;
}

void ProbSetCheck::initialise( std::size_t nSamples, std::size_t nAlleles ) {
	REQUIRE( m_state == eNone ) ;
	REQUIRE( nSamples == m_number_of_samples ) ;
	REQUIRE( nAlleles == 2 ) ;
	m_state = eSetNumberOfSamples ;
}

bool ProbSetCheck::set_sample( std::size_t i ) {
	REQUIRE( i < m_number_of_samples ) ;
	bool const correctState = (m_state == eSetNumberOfSamples) || (m_state == eSetValue) ;
	REQUIRE( correctState ) ;
	REQUIRE( m_entry_i == m_number_of_entries ) ;
	m_sample_i = i ;
	m_state = eSetSample ;
	return true ;
}

void ProbSetCheck::set_number_of_entries( std::size_t n, genfile::OrderType const order_type, genfile::ValueType const value_type ) {
#if DEBUG > 2
	std::cerr << "ProbSetApiCheck::set_number_of_entries(): n = " << n
		<< ", order_type = " << order_type << ".\n" ;
#endif
	REQUIRE( m_state == eSetSample ) ;
	bool const correctOrderType = (
		( order_type == genfile::ePerUnorderedGenotype && n == 3)
		||
		( order_type == genfile::ePerPhasedHaplotypePerAllele && n == 4)
	) ;
	REQUIRE( correctOrderType ) ;
	REQUIRE( value_type == genfile::eProbability ) ;
	m_number_of_entries = n ;
	m_order_type = order_type ;
	m_entry_i = 0 ;
	m_state = eSetNumberOfEntries ;
}

void ProbSetCheck::operator()( genfile::MissingValue const value ) {
	REQUIRE(0) ;
}

void ProbSetCheck::operator()( double const value ) {
	bool const correctEntry = (
		(m_entry_i == 0 && m_state == eSetNumberOfEntries)
		 || (m_entry_i > 0 && m_state == eSetValue)
	) ;
	REQUIRE( correctEntry ) ;
	REQUIRE( m_entry_i < m_number_of_entries ) ;
#if DEBUG > 2
	std::cerr << format( "ProbabilitySetter: sample %d, entry %d of %d.\n", m_sample_i, m_entry_i, m_number_of_entries ) ;
#endif
	//REQUIRE( value == Approx( m_get_expected_probs( m_sample_i, m_entry_i ) ) ;
	++m_entry_i ;
	m_state = eSetValue ;
}
