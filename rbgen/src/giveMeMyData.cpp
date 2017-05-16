#include "genfile/bgen/View.hpp"
#include "genfile/bgen/IndexQuery.hpp"
#include <Rcpp.h>

#define DEBUG 0

namespace {
	struct DataFiller {
		typedef Rcpp::IntegerVector IntegerVector ;
		typedef Rcpp::NumericVector NumericVector ;
		typedef Rcpp::Dimension Dimension ;
		
		DataFiller(
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
		{}
		
		// Called once allowing us to set storage.
		void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
			// Nothing to do, but check storage is big enough.
			if(!(
				(number_of_alleles == 2)
				&& ( m_data_dimension[0] > m_variant_i )
				&& ( m_data_dimension[1] == number_of_samples )
				&& ( m_data_dimension[2] == 3 )
				&& ( m_ploidy_dimension[0] > m_variant_i )
				&& ( m_ploidy_dimension[1] == number_of_samples )
			)) {
#if DEBUG
				std::cerr << "initialise(): invalid constraint.\n" ;
				std::cerr << "Data dimension:\n" << m_data_dimension[0] << ", " << m_data_dimension[1] << ", " << m_data_dimension[2] << "\n" ;
				std::cerr << "Ploidy dimension:\n" << m_ploidy_dimension[0] << ", " << m_ploidy_dimension[1] << "\n" ;
				std::cerr << "Number of samples: " << number_of_samples << "\n" ;
				std::cerr << "Number of alleles: " << number_of_alleles << "\n" ;
				std::cerr << "Variant: " << m_variant_i << "\n" ;
#endif
				throw genfile::bgen::BGenError() ;
			}
		}
	
		// If present with this signature, called once after initialise()
		// to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
		// This enables us to set up storage for the data ahead of time.
		void set_min_max_ploidy( uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries ) {
			assert( min_ploidy == 2 ) ;
			assert( max_ploidy == 2 ) ;
			// Now to do here given our assumptions
			m_sample_i = 0 ;
		}
	
		// Called once per sample to determine whether we want data for this sample
		bool set_sample( std::size_t i ) {
			m_sample_i = i ;
			// Yes, here we want info for all samples.
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
				throw genfile::bgen::BGenError() ;
			}
			// assert( value_type == genfile::eProbability ) ;
			m_order_type = order_type ;
			int flatIndex = m_variant_i + m_sample_i * m_ploidy_dimension[0] ;
			(*m_ploidy)[ flatIndex ] = ploidy ;
		}

		void set_value( uint32_t entry_i, double value ) {
			int flatIndex = m_variant_i + m_sample_i * m_data_dimension[0] + entry_i * m_data_dimension[0] * m_data_dimension[1] ;
#if DEBUG
			std::cerr << "Setting data for index " << m_variant_i << ", " << m_sample_i << ", " << entry_i << ": index " << flatIndex << "...\n" << std::flush ;
#endif
			(*m_data)[ flatIndex ] = value ;
		}

		void set_value( uint32_t entry_i, genfile::MissingValue value ) {
			int flatIndex = m_variant_i + m_sample_i * m_data_dimension[0] + entry_i * m_data_dimension[0] * m_data_dimension[1] ;
#if DEBUG
			std::cerr << "Setting data for index " << m_variant_i << ", " << m_sample_i << ", " << entry_i << ": index " << flatIndex << "...\n" << std::flush ;
#endif
			(*m_data)[ flatIndex ] = NA_REAL ;
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
Rcpp::List giveMeMyData(
	std::string const& filename,
	std::string const& chromosome,
	uint32_t start,
	uint32_t end
) {
	using namespace genfile::bgen ;
	using namespace Rcpp ;

	View::UniquePtr view = View::create( filename ) ;
	
	IndexQuery::UniquePtr query = IndexQuery::create( filename	+ ".bgi" ) ;
	{
	query->include_range( IndexQuery::GenomicRange( chromosome, start, end )) ;
	query->initialise() ;
	view->set_query( query ) ;
	}

	List result ;
	{
		std::size_t const NV = view->number_of_variants() ;
		CharacterVector chromosomes( NV ) ;
		IntegerVector positions( NV ) ;
		CharacterVector rsids( NV ) ;
		CharacterVector first_alleles( NV ) ;
		CharacterVector alternate_alleles( NV ) ;
	
		std::string SNPID ;
		std::string rsid ;
		std::string chromosome ;
		uint32_t position ;
		std::vector< std::string > alleles ;
	
		Dimension ploidy_dimension = Dimension( NV, view->number_of_samples() ) ;
		IntegerVector ploidy = IntegerVector( ploidy_dimension ) ;
		Dimension data_dimension = Dimension( NV, view->number_of_samples(), 3ul ) ;
		NumericVector data = NumericVector( data_dimension ) ;
		
		StringVector sampleNames( view->number_of_samples() ) ;
		view->get_sample_ids( set_sample_names( &sampleNames )) ;
		
#if DEBUG
		std::cerr << "ploidy.size() == " << ploidy.size() << "!\n" ;
		std::cerr << "data.size() == " << data.size() << "!\n" ;
#endif
		
		for( std::size_t i = 0 ; i < NV; ++i ) {
			DataFiller filler( &ploidy, ploidy_dimension, &data, data_dimension, i ) ;
			view->read_variant( &SNPID, &rsid, &chromosome, &position, &alleles ) ;
			view->read_genotype_data_block( filler ) ;
			//view->ignore_genotype_data_block() ;
			rsids[i] = rsid ;
			first_alleles[i] = alleles[0] ;
			alternate_alleles[i] = alleles[1] ;
			chromosomes[i] = chromosome ;
			positions[i] = position ;
		}

		data[0] = 1 ;
		data[1] = 2 ;
		data[3] = 3 ;

		// For this example we assume diploid samples and two alleles
		DataFrame variants = DataFrame::create(
			Named("chromosome") = chromosomes,
			Named("position") = positions,
			Named("rsid") = rsids,
			Named("allele0") = first_alleles,
			Named("allele1") = alternate_alleles
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
	
		data.attr( "dimnames" ) = dataNames ;
	
		result[ "variants" ] = variants ;
		result[ "samples" ] = sampleNames ;
		result[ "ploidy" ] = ploidy ;
		result[ "data" ] = data ;
	}
	return( result ) ;
}
