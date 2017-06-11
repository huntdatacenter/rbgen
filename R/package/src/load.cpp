#include "genfile/bgen/View.hpp"
#include "genfile/bgen/IndexQuery.hpp"
#include <sstream>
#include <Rcpp.h>

// #define DEBUG 1

namespace {
	template< typename Integer >
	std::string atoi( Integer const value ) {
		std::ostringstream stream ;
		stream << value ;
		return stream.str() ;
	}
	
	struct DataSetter {
		typedef Rcpp::IntegerVector IntegerVector ;
		typedef Rcpp::NumericVector NumericVector ;
		typedef Rcpp::Dimension Dimension ;

		DataSetter(
			IntegerVector* ploidy,
			Dimension const& ploidy_dimension,
			NumericVector* data,
			Dimension const& data_dimension,
			std::size_t variant_i
		):
			m_ploidy( ploidy ),
			m_ploidy_dimension( ploidy_dimension ),
			m_data( data ),
			m_data_dimension( data_dimension ),
			m_variant_i( variant_i )
		{
			assert( m_data_dimension[0] == m_ploidy_dimension[0] ) ;
			assert( m_data_dimension[1] == m_ploidy_dimension[1] ) ;
			assert( m_variant_i < m_data_dimension[0] ) ;
			assert( m_data_dimension[2] >= 3 ) ;
		}
	
		// Called once allowing us to set storage.
		void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
			// Nothing to do, but check storage is big enough.
			if( number_of_alleles != 2 ) {
				throw std::invalid_argument( "number_of_samples != 2" ) ;
			}
			if( m_data_dimension[1] != number_of_samples ) {
				throw std::invalid_argument( "Wrong number of samples ("
					+ atoi( number_of_samples ) + ", expected " + atoi( m_data_dimension[1] ) + ")" ) ;
			}
		}

		// If present with this signature, called once after initialise()
		// to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
		// This enables us to set up storage for the data ahead of time.
		void set_min_max_ploidy( genfile::bgen::uint32_t min_ploidy, genfile::bgen::uint32_t max_ploidy, genfile::bgen::uint32_t min_entries, genfile::bgen::uint32_t max_entries ) {
			if( max_ploidy > 2 ) {
				throw std::invalid_argument( "max_ploidy=" + atoi( max_ploidy ) + " (expected at most 2)" ) ;
			}
		}

		// Called once per sample to determine whether we want data for this sample
		bool set_sample( std::size_t i ) {
			m_sample_i = i ;
			return true ;
		}

		// Called once per sample to set the number of probabilities that are present.
		void set_number_of_entries(
			std::size_t ploidy,
			std::size_t number_of_entries,
			genfile::OrderType order_type,
			genfile::ValueType value_type
		) {
			if( value_type != genfile::eProbability ) {
				throw std::invalid_argument(
					"value_type ("
					+ atoi( value_type ) + ", expected "
					+ atoi( genfile::eProbability ) + "=genfile::eProbability)"
				) ;
			}
			m_order_type = order_type ;
			int const index = m_variant_i + m_sample_i * m_ploidy_dimension[0] ;
			(*m_ploidy)[ index ] = ploidy ;
		}

		void set_value( genfile::bgen::uint32_t entry_i, double value ) {
			int const index = m_variant_i + m_sample_i * m_data_dimension[0] + entry_i * m_data_dimension[0] * m_data_dimension[1] ;
#if DEBUG
			std::cerr << "Setting data for index " << m_variant_i << ", " << m_sample_i << ", " << entry_i << ": index " << index << "...\n" << std::flush ;
#endif
			(*m_data)[ index ] = value ;
		}

		void set_value( genfile::bgen::uint32_t entry_i, genfile::MissingValue value ) {
			int const index = m_variant_i + m_sample_i * m_data_dimension[0] + entry_i * m_data_dimension[0] * m_data_dimension[1] ;
#if DEBUG
			std::cerr << "Setting data for index " << m_variant_i << ", " << m_sample_i << ", " << entry_i << ": index " << index << "...\n" << std::flush ;
#endif
			(*m_data)[ index ] = NA_REAL ;
		}

		void finalise() {
		// nothing to do
		}

	private:
		Rcpp::IntegerVector* m_ploidy ;
		Rcpp::Dimension const m_ploidy_dimension ;
		Rcpp::NumericVector* m_data ;
		Rcpp::Dimension const m_data_dimension ;

		std::size_t const m_variant_i ;
		std::size_t m_sample_i ;

		genfile::OrderType m_order_type ;
	} ;

	struct set_sample_names {
		set_sample_names( Rcpp::StringVector* result ):
			m_result( result ),
			m_index(0)
		{
			assert( result != 0 ) ;
		}
		
		void operator()( std::string const& value ) {
			(*m_result)[m_index++] = value ;
		}
	private:
		Rcpp::StringVector* m_result ;
		std::size_t m_index ;
	} ;
}

// [[Rcpp::export]]
Rcpp::List load(
	std::string const& filename,
	Rcpp::DataFrame const& ranges
) {
	using namespace genfile::bgen ;
	using namespace Rcpp ;

	View::UniquePtr view = View::create( filename ) ;
	{
		StringVector const& chromosome = ranges["chromosome"] ;
		IntegerVector const& start = ranges["start"] ;
		IntegerVector const& end = ranges["end"] ;
		
		IndexQuery::UniquePtr query = IndexQuery::create( filename + ".bgi" ) ;
		for( int i = 0; i < ranges.nrows(); ++i ) {
			query->include_range( IndexQuery::GenomicRange( std::string( chromosome[i] ), start[i], end[i] )) ;
		}
		query->initialise() ;
		view->set_query( query ) ;
	}

	std::size_t const number_of_variants = view->number_of_variants() ;
	std::size_t const number_of_samples = view->number_of_samples() ;

	// For this example we assume diploid samples and two alleles
	StringVector chromosomes( number_of_variants ) ;
	IntegerVector positions( number_of_variants ) ;
	StringVector rsids( number_of_variants ) ;
	StringVector allele0s( number_of_variants ) ;
	StringVector allele1s( number_of_variants ) ;
	StringVector sampleNames( number_of_samples ) ;

	Dimension data_dimension = Dimension( number_of_variants, number_of_samples, 3ul ) ;
	Dimension ploidy_dimension = Dimension( number_of_variants, number_of_samples ) ;

	NumericVector data = NumericVector( data_dimension ) ;
	IntegerVector ploidy = IntegerVector( ploidy_dimension ) ;

	view->get_sample_ids( set_sample_names( &sampleNames ) ) ;

	std::string SNPID, rsid, chromosome ;
	genfile::bgen::uint32_t position ;
	std::vector< std::string > alleles ;

	for( std::size_t variant_i = 0; variant_i < number_of_variants; ++variant_i ) {
		view->read_variant( &SNPID, &rsid, &chromosome, &position, &alleles ) ;
		chromosomes[variant_i] = chromosome ;
		positions[variant_i] = position ;
		rsids[variant_i] = rsid ;
		allele0s[variant_i] = alleles[0] ;
		allele1s[variant_i] = alleles[1] ;

	DataSetter setter(
		&ploidy, ploidy_dimension,
		&data, data_dimension,
		variant_i
	) ;
		view->read_genotype_data_block( setter ) ; // will be fixed later
	}

	DataFrame variants = DataFrame::create(
		Named("chromosome") = chromosomes,
		Named("position") = positions,
		Named("rsid") = rsids,
		Named("allele0") = allele0s,
		Named("allele1") = allele1s
	) ;
	variants.attr( "row.names" ) = rsids ;

	StringVector genotypeNames(3) ;
	genotypeNames[0] = "0" ;
	genotypeNames[1] = "1" ;
	genotypeNames[2] = "2" ;

	List dataNames = List(3) ;
	dataNames[0] = rsids ;
	dataNames[1] = sampleNames ;
	dataNames[2] = genotypeNames ;

	List ploidyNames = List(2) ;
	ploidyNames[0] = rsids ;
	ploidyNames[1] = sampleNames ;

	data.attr( "dimnames" ) = dataNames ;
	ploidy.attr( "dimnames" ) = ploidyNames ;

	List result ;
	result[ "variants" ] = variants ;
	result[ "samples" ] = sampleNames ;
	result[ "ploidy" ] = ploidy ;
	result[ "data" ] = data ;

	return( result ) ;
}
