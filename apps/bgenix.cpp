
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <fstream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>
#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "genfile/bgen/bgen.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "../../bgen_revision_autogenerated.hpp"

namespace bfs = boost::filesystem ;

// #define DEBUG 1

namespace globals {
	std::string const program_name = "bgenix" ;
	std::string const program_version = bgen_revision ;
}

struct IndexBgenOptionProcessor: public appcontext::CmdLineOptionProcessor
{
public:
	std::string get_program_name() const { return globals::program_name ; }

	void declare_options( appcontext::OptionProcessor& options ) {
		// Meta-options
		options.set_help_option( "-help" ) ;

		options.declare_group( "Input / output file options" ) ;
		options[ "-g" ]
			.set_description(
				"Path of bgen file to operate on.  (An optional form where \"-g\" is omitted and the filename is specified as the first argument, i.e. bgenix <filename>, can also be used)."
			)
			.set_takes_single_value()
			.set_takes_value_by_position(1)
			.set_is_required()
		;

		options[ "-clobber" ]
			.set_description(
				"Specify that bgenix should overwrite existing index file if it exists."
			)
		;
		
		options[ "-table" ]
			.set_description( "Specify the table (or view) that bgenix should read the file index from."
				"This only affects reading the index file.  The named table or view should have the"
				"same schema as the Variant table written by bgenix on index creation."
			)
			.set_takes_single_value()
			.set_default_value( "Variant" )
		;
		
		options.declare_group( "Variant selection options" ) ;
		options[ "-incl-range" ]
			.set_description(
				"Include variants in the specified genomic interval in the output. "
				" The argument must be of the form <chr>:<pos1>-<pos2> where <chr> is a chromosome identifier "
				" and pos1 and pos2 are positions with pos2 >= pos1. "
				" At most one of pos1 and pos2 can also be omitted, in which case the range extends to the start or"
				" end of the chromosome as appropriate. "
				" Position ranges are treated as closed (i.e. <pos1> and <pos2> are included in the range)."
			)
			.set_takes_value_by_position(2)
			.set_takes_values_until_next_option()
		;

		options[ "-excl-range" ]
			.set_description(
				"Exclude variants in the specified genomic interval from the output. "
				"See the description of -incl-range for details."
				"If this is specified multiple times, variants in any of the specified ranges will be excluded."
			)
			.set_takes_values_until_next_option()
		;

		options[ "-incl-rsids" ]
			.set_description(
				"Include variants with the specified rsid(s) in the output. "
			)
			.set_takes_values_until_next_option()
		;

		options[ "-excl-rsids" ]
			.set_description(
				"Exclude variants with the specified rsid(s) from the output. "
			)
			.set_takes_values_until_next_option()
		;
		
		options.declare_group( "Output options" ) ;
		options[ "-list" ]
			.set_description( "Suppress BGEN output; instead output a list of variants." ) ;
		options[ "-with-rowid" ]
			.set_description( "Create an index file without using the 'WITHOUT ROWID' tables.  These are suitable for use with sqlite versions < 3.8.2" ) ;
	}
} ;

struct BgenProcessor {
	
	BgenProcessor( std::string const& filename ):
		m_filename( filename ),
		m_state( e_NotOpen )
	{
		// Open the stream
		m_stream.reset(
			new std::ifstream( filename, std::ifstream::binary )
		) ;
		if( !*m_stream ) {
			throw std::invalid_argument( filename ) ;
		}
		m_state = e_Open ;

		// Read the offset, header, and sample IDs if present.
		genfile::bgen::read_offset( *m_stream, &m_offset ) ;
		genfile::bgen::read_header_block( *m_stream, &m_context ) ;

		// Skip anything else until the first variant
		m_stream->seekg( m_offset + 4 ) ;

		// We update track of the start and end of each variant
		m_current_variant_position = ( m_offset + 4 ) ;

		// We keep track of state (though it's not really needed for this implementation.)
		m_state = e_ReadyForVariant ;
	}

	std::streampos current_position() const {
		return m_current_variant_position ;
	}
	
	void seek_to( std::streampos pos ) {
		m_current_variant_position = pos ;
		m_stream->seekg( pos ) ;
		m_state = e_ReadyForVariant ;
	}

	uint32_t number_of_variants() const {
		return m_context.number_of_variants ;
	}

	// Attempt to read identifying information about a variant from the bgen file, returning
	// it in the given fields.
	// If this method returns true, data was successfully read, and it should be safe to call read_probs()
	// or ignore_probs().
	// If this method returns false, data was not successfully read indicating the end of the file.
	bool read_variant(
		std::string* SNPID,
		std::string* rsid,
		std::string* chromosome,
		uint32_t* position,
		std::vector< std::string >* alleles
	) {
		assert( m_state == e_ReadyForVariant ) ;
		
		if(
			genfile::bgen::read_snp_identifying_data(
				*m_stream, m_context,
				SNPID, rsid, chromosome, position,
				[&alleles]( std::size_t n ) { alleles->resize( n ) ; },
				[&alleles]( std::size_t i, std::string const& allele ) { alleles->at(i) = allele ; }
			)
		) {
			m_state = e_ReadyForProbs ;
			return true ;
		} else {
			return false ;
		}
	}
	
	// Ignore genotype probability data for the SNP just read using read_variant()
	// After calling this method it should be safe to call read_variant()
	// to fetch the next variant from the file.
	void ignore_probs() {
		genfile::bgen::ignore_genotype_data_block( *m_stream, m_context ) ;
		m_current_variant_position = m_stream->tellg() ;
		m_state = e_ReadyForVariant ;
	}

private:
	std::string const m_filename ;
	std::unique_ptr< std::istream > m_stream ;

	// bgen::Context object holds information from the header block,
	// including bgen flags
	genfile::bgen::Context m_context ;

	// offset byte from top of bgen file.
	uint32_t m_offset ;

	// We keep track of our state in the file.
	// Not strictly necessary for this implentation but makes it clear that
	// calls must be read_variant() followed by read_probs() (or ignore_probs())
	// repeatedly.
	enum State { e_NotOpen = 0, e_Open = 1, e_ReadyForVariant = 2, e_ReadyForProbs = 3, eComplete = 4 } ;
	State m_state ;
	
	// To avoid issues with tellg() and failbit, we store the stream position at suitable points
	std::streampos m_current_variant_position ;
} ;

struct IndexBgenApplication: public appcontext::ApplicationContext
{
public:
	IndexBgenApplication( int argc, char** argv ):
		appcontext::ApplicationContext(
			globals::program_name,
			globals::program_version,
			std::auto_ptr< appcontext::OptionProcessor >( new IndexBgenOptionProcessor ),
			argc,
			argv,
			"-log"
		)
	{
		setup() ;
	}

private:
	std::string m_bgen_filename ;
	std::string m_index_filename ;

private:
	void setup() {
		m_bgen_filename = options().get< std::string >( "-g" ) ;
		m_index_filename = m_bgen_filename + ".bgi" ;
		if( !bfs::exists( m_bgen_filename )) {
			throw std::invalid_argument( m_bgen_filename ) ;
		}
		if( options().check( "-clobber" ) || !bfs::exists( m_index_filename )) {
			create_bgen_index( m_bgen_filename, m_index_filename ) ;
		}
		if( options().check_if_option_was_supplied_in_group( "Variant selection options" )) {
			process_selection( m_bgen_filename, m_index_filename ) ;
		}
	}


	void create_bgen_index( std::string const& bgen_filename, std::string const& index_filename ) {
		db::Connection::UniquePtr result ;
		ui().logger()
			<< boost::format( "%s: creating index for \"%s\" in \"%s\"...\n" ) % globals::program_name % bgen_filename % index_filename ;

		try {
			assert( !bfs::exists( index_filename ) ) ;
			assert( !bfs::exists( index_filename + ".tmp" ) ) ;
			result = create_bgen_index_unsafe( bgen_filename, index_filename + ".tmp" ) ;
			bfs::rename( index_filename + ".tmp", index_filename ) ;
		} catch( db::StatementStepError const& e ) {
			ui().logger() << "!! Error in \"" << e.spec() << "\": " << e.description() << ".\n" ;
			bfs::remove( index_filename + ".tmp" ) ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		} catch( ... ) {
			// Remove the incomplete attempt at an index file.
			bfs::remove( index_filename + ".tmp" ) ;
			throw ;
		}
	}
	
	db::Connection::UniquePtr create_bgen_index_unsafe( std::string const& bgen_filename, std::string const& index_filename ) {
		db::Connection::UniquePtr connection = db::Connection::create( index_filename ) ;

		connection->run_statement( "PRAGMA locking_mode = EXCLUSIVE ;" ) ;
		connection->run_statement( "PRAGMA journal_mode = MEMORY ;" ) ;
		
		db::Connection::ScopedTransactionPtr transaction = connection->open_transaction( 240 ) ;
		setup_index_file( *connection ) ;
		// Close and open the transaction
		transaction.reset() ;

		db::Connection::StatementPtr insert_variant_stmt = connection->get_statement(
			"INSERT INTO Variant( chromosome, position, rsid, number_of_alleles, allele1, allele2, file_start_position, size_in_bytes ) "
			"VALUES( ?, ?, ?, ?, ?, ?, ?, ? )"
		) ;

		BgenProcessor processor( bgen_filename ) ;

		ui().logger()
			<< boost::format( "%s: Opened \"%s\" with %d variants...\n" ) % globals::program_name % bgen_filename % processor.number_of_variants() ;
		
		std::string chromosome, rsid, SNPID ;
		uint32_t position ;
		std::vector< std::string > alleles ;
		alleles.reserve(100) ;
		std::size_t const chunk_size = 10 ;
		
		transaction = connection->open_transaction( 240 ) ;
		
		{
			auto progress_context = ui().get_progress_context( "Building BGEN index" ) ;
			std::size_t variant_count = 0;
			int64_t file_pos = int64_t( processor.current_position() ) ;
			try {
				while( processor.read_variant( &SNPID, &rsid, &chromosome, &position, &alleles ) ) {
					processor.ignore_probs() ;
					int64_t file_end_pos = int64_t( processor.current_position() ) ;
					assert( alleles.size() > 1 ) ;
					assert( (file_end_pos - file_pos) > 0 ) ;
#if DEBUG
					std::cerr << chromosome << " " << position << " " << rsid << " " << file_pos << ".\n" ;
#endif
					insert_variant_stmt
						->bind( 1, chromosome )
						.bind( 2, position )
						.bind( 3, rsid )
						.bind( 4, int64_t( alleles.size() ) )
						.bind( 5, alleles[0] )
						.bind( 6, alleles[1] )
						.bind( 7, file_pos )
						.bind( 8, file_end_pos - file_pos )
						.step()
					;
					insert_variant_stmt->reset() ;
				
					progress_context( ++variant_count, processor.number_of_variants() ) ;
			
				// Make sure and commit every 10000 SNPs.
					if( variant_count % chunk_size == 0 ) {
		//				ui().logger()
		//					<< boost::format( "%s: Writing variants %d-%d...\n" ) % processor.number_of_variants() % (variant_count-chunk_size) % (variant_count-1) ;
					
						transaction.reset() ;
						transaction = connection->open_transaction( 240 ) ;
					}
					file_pos = file_end_pos ;
				}
			}
			catch( db::StatementStepError const& e ) {
				ui().logger() << "Last observed variant was " << SNPID << " " << rsid << " " << chromosome << " " << position ;
				for( std::size_t i = 0; i < alleles.size(); ++i ) {
					ui().logger() << " " << alleles[i] ;
				}
				ui().logger() << "\n" ;
				throw ;
			}
		}
		return connection ;
	}
	
	void setup_index_file( db::Connection& connection ) {
		std::string const tag = options().check( "-with-rowid" ) ? "" : " WITHOUT ROWID" ;

		connection.run_statement(
			"CREATE TABLE Variant ("
			"  chromosome TEXT NOT NULL,"
			"  position INT NOT NULL,"
			"  rsid TEXT NOT NULL,"
			"  number_of_alleles INT NOT NULL,"
			"  allele1 TEXT NOT NULL,"
			"  allele2 TEXT NULL,"
			"  file_start_position INT NOT NULL," // 
			"  size_in_bytes INT NOT NULL,"       // We put these first to minimise cost of retrieval
			"  PRIMARY KEY (chromosome, position, rsid, allele1, allele2, file_start_position )"
			")" + tag
		) ;
	}
	
	void process_selection( std::string const& bgen_filename, std::string const& index_filename ) const {
		try {
			process_selection_unsafe( bgen_filename, index_filename ) ;
		}
		catch( ... ) {
			throw ;
		}
	}

	void process_selection_unsafe( std::string const& bgen_filename, std::string const& index_filename ) const {
		db::Connection::UniquePtr connection = db::Connection::create( index_filename, "read" ) ;
		std::string const& select_sql = get_select_sql( *connection ) ;
#if DEBUG
		std::cerr << "Running sql: \"" << select_sql << "\"...\n" ;
#endif
		std::vector< std::pair< int64_t, int64_t> > positions ;
		{
			auto progress_context = ui().get_progress_context( "Selecting variants" ) ;
			progress_context( 0, boost::optional< std::size_t >() ) ;
			db::Connection::StatementPtr stmt = connection->get_statement( select_sql ) ;
			{
				positions.reserve( 1000000 ) ;
				std::size_t batch_i = 0 ;
				for( stmt->step() ; !stmt->empty(); stmt->step(), ++batch_i ) {
					int64_t const pos = stmt->get< int64_t >( 0 ) ;
					int64_t const size = stmt->get< int64_t >( 1 ) ;
					assert( pos >= 0 ) ;
					assert( size >= 0 ) ;
					positions.push_back( std::make_pair( int64_t( pos ), int64_t( size ))) ;
					progress_context( positions.size(), boost::optional< std::size_t >() ) ;
				}
			}
		}

		if( options().check( "-list" ) ) {
			process_selection_list( bgen_filename, positions ) ;
		} else {
			process_selection_bgen( bgen_filename, positions ) ;
		}
	}

	void process_selection_bgen( std::string const& bgen_filename, std::vector< std::pair< int64_t, int64_t> > const& positions ) const {
		std::ifstream bgen_file( bgen_filename, std::ios::binary ) ;
		uint32_t offset = 0 ;

		using namespace genfile ;
		bgen::Context context ;
		bgen::read_offset( bgen_file, &offset ) ;
		bgen::read_header_block( bgen_file, &context ) ;

		// Write the new context after adjusting the variant count.
		context.number_of_variants = positions.size() ;
		bgen::write_offset( std::cout, offset ) ;
		bgen::write_header_block( std::cout, context ) ;
		std::istreambuf_iterator< char > inIt( bgen_file ) ;
		std::istreambuf_iterator< char > endInIt ;
		std::ostreambuf_iterator< char > outIt( std::cout ) ;

		// Copy everything else up to the start of the data
		std::copy_n( inIt, offset - context.header_size(), outIt ) ;

		{
			auto progress_context = ui().get_progress_context( "Processing " + std::to_string( positions.size() ) + " variants" ) ;
			// Now we go for it
			for( std::size_t i = 0; i < positions.size(); ++i ) {
				bgen_file.seekg( positions[i].first ) ;
				std::istreambuf_iterator< char > inIt( bgen_file ) ;
				std::copy_n( inIt, positions[i].second, outIt ) ;
				progress_context( i+1, positions.size() ) ;
			}
		}
		std::cerr << boost::format( "%s: wrote data for %d variants to stdout.\n" ) % globals::program_name % positions.size() ;
	}
	
	void process_selection_list( std::string const& bgen_filename, std::vector< std::pair< int64_t, int64_t> > const& positions ) const {
		BgenProcessor processor( bgen_filename ) ;
		std::cout << boost::format( "# %s: started %s\n" ) % globals::program_name % appcontext::get_current_time_as_string() ;
		std::cout << "alternate_ids\trsid\tchromosome\tposition\tnumber_of_alleles\tfirst_allele\talternative_alleles\n" ;
		
		std::string SNPID, rsid, chromosome ;
		uint32_t position ;
		std::vector< std::string > alleles ;

		for( std::size_t i = 0; i < positions.size(); ++i ) {
			processor.seek_to( positions[i].first ) ;
			bool success = processor.read_variant(
				&SNPID, &rsid, &chromosome, &position, &alleles
			) ;
			std::cout << SNPID << "\t" << rsid << "\t" << chromosome << "\t" << position << "\t" << alleles.size() << "\t" << alleles[0] << "\t" ;
			for( std::size_t allele_i = 1; allele_i < alleles.size(); ++allele_i ) {
				std::cout << (( allele_i > 1 ) ? "," : "" ) << alleles[allele_i] ;
			}
			std::cout << "\n" ;
			if( !success ) {
				throw std::invalid_argument( "positions" ) ;
			}
			processor.ignore_probs() ;
		}
		
		std::cout << boost::format( "# %s: success, total %d variants.\n" ) % globals::program_name % positions.size() ;
	}

	boost::tuple< std::string, uint32_t, uint32_t > parse_range( std::string const& spec ) const {
		std::size_t colon_pos = spec.find( ':' ) ;
		if ( colon_pos == std::string::npos ) {
			throw std::invalid_argument( "spec=\"" + spec + "\"" ) ;
		}

		std::vector< std::string > pieces ;
		pieces.push_back( spec.substr( 0, colon_pos )) ;
		pieces.push_back( spec.substr( colon_pos+1, spec.size() )) ;

		if ( pieces.size() != 2 ) {
			throw std::invalid_argument( "spec=\"" + spec + "\"" ) ;
		}

		std::size_t separator_pos = pieces[1].find( '-' ) ;
		if ( separator_pos == std::string::npos ) {
			throw std::invalid_argument( "spec=\"" + spec + "\"" ) ;
		}

		std::string chromosome( pieces[0] ) ;
		int pos1 = (separator_pos == 0) ? 0 : std::stoi( pieces[1].substr( 0, separator_pos ) ) ;
		int pos2 = (separator_pos == (pieces[1].size()-1)) ? std::numeric_limits< int >::max() : std::stoi( pieces[1].substr( separator_pos + 1, pieces[1].size() ) ) ;
		assert( pos1 >= 0 ) ;
		assert( pos2 >= pos1 ) ;

		return boost::make_tuple( chromosome, uint32_t( pos1 ), uint32_t( pos2 ) ) ;
	}

	std::string get_select_sql( db::Connection& connection ) const {
		std::string result = "SELECT file_start_position, size_in_bytes FROM `"
			+ options().get_value("-table")
			+ "` V" ;
		std::string join ;
		std::string inclusion = "" ;
		std::string exclusion = "" ;
		if( options().check( "-incl-range" )) {
			auto elts = options().get_values< std::string >( "-incl-range" ) ;
			for( std::string const& elt: elts ) {
				boost::tuple< std::string, uint32_t, uint32_t > range = parse_range( elt ) ;
				inclusion += ((inclusion.size() > 0) ? " OR " : "" ) + (
					boost::format( "( chromosome == '%s' AND position BETWEEN %d AND %d )" ) % range.get<0>() % range.get<1>() % range.get<2>()
				).str() ;
			}
		}
		if( options().check( "-excl-range" )) {
			auto elts = options().get_values< std::string >( "-excl-range" ) ;
			for( std::string const& elt: elts ) {
				boost::tuple< std::string, uint32_t, uint32_t > range = parse_range( elt ) ;
				exclusion += ( exclusion.size() > 0 ? " AND" : "" ) + (
					boost::format( " NOT ( chromosome == '%s' AND position BETWEEN %d AND %d )" ) % range.get<0>() % range.get<1>() % range.get<2>()
				).str() ;
			}
		}
		if( options().check( "-incl-rsids" )) {
			auto const ids = collect_unique_ids( options().get_values< std::string >( "-incl-rsids" ));
			connection.run_statement( "CREATE TEMP TABLE tmpIncludedId( rsid TEXT NOT NULL PRIMARY KEY ) WITHOUT ROWID" ) ;
			db::Connection::StatementPtr insert_stmt = connection.get_statement( "INSERT INTO tmpIncludedId( rsid ) VALUES( ? )" ) ;
			for( auto elt: ids ) {
				insert_stmt
					->bind( 1, elt )
					.step() ;
				insert_stmt->reset() ;
			}
			join += " INNER JOIN tmpIncludedId T1 ON T1.rsid == V.rsid" ;
		}
		if( options().check( "-excl-rsids" )) {
			auto const ids = collect_unique_ids( options().get_values< std::string >( "-incl-rsids" ));
			connection.run_statement( "CREATE TEMP TABLE tmpExcludedId( rsid TEXT NOT NULL PRIMARY KEY ) WITHOUT ROWID" ) ;
			db::Connection::StatementPtr insert_stmt = connection.get_statement( "INSERT INTO tmpExcludedId( rsid ) VALUES( ? )" ) ;
			for( auto elt: ids ) {
				insert_stmt
					->bind( 1, elt )
					.step() ;
				insert_stmt->reset() ;
			}
			join += " LEFT OUTER JOIN tmpExcludedId TE ON TE.rsid == V.rsid" ;
			exclusion += ( exclusion.size() > 0 ? " AND" : "" ) + std::string( " TE.rsid IS NULL" ) ;
		}
		inclusion = "(" + inclusion + ")" ;
		exclusion += "(" + exclusion + ")" ;
		std::string const where = "WHERE " + inclusion + " AND " + exclusion ;
		std::string const orderBy = "ORDER BY chromosome, position, rsid, allele1, allele2" ;
		return result + " " + join + " " + where +  inclusion + " AND " + exclusion + " " + orderBy ;
	}
	
	std::vector< std::string > collect_unique_ids( std::vector< std::string > const& ids_or_filenames ) const {
		std::vector< std::string > result ;
		for( auto elt: ids_or_filenames ) {
			if( bfs::exists( elt )) {
				std::ifstream f( elt ) ;
				std::copy(
					std::istream_iterator< std::string >( f ),
					std::istream_iterator< std::string >(),
					std::back_inserter< std::vector< std::string > >( result )
				) ;
			} else {
				result.push_back( elt ) ;
			}
		}
		// now sort and uniqueify them...
		std::sort( result.begin(), result.end() ) ;
		std::vector< std::string >::const_iterator newBack = std::unique( result.begin(), result.end() ) ;
		result.resize( newBack - result.begin() ) ;
		return result ;
	}

	
} ;

int main( int argc, char** argv ) {
    try {
		IndexBgenApplication app( argc, argv ) ;
    }
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}
