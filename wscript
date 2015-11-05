import glob

srcdir="."
APPNAME = "bgen"
VERSION = "1.2"

def options( opt ):
	opt.load( 'compiler_cxx' )
	
def configure( cfg ):
	cfg.load( 'compiler_cxx')
	cfg.env.CXXFLAGS = [ '-g' ]
#	cfg.env.CFLAGS = [ '-O3' ]
	cfg.check_cxx( lib='z', uselib_store='zlib', msg = 'zlib' )

def build( bld ):
	bld.stlib(
		source = bld.path.ant_glob( 'src/*.cpp' ),
		target = 'bgen',
		includes = 'genfile/include',
		export_includes = 'genfile/include'
	)

	bld.program(
		source = 'example/bgen_to_vcf.cpp',
		target = 'bgen_to_vcf',
		use = 'bgen zlib',
		cxxflags = ['-std=c++11']
	)

	bld.program(
		source = bld.path.ant_glob( 'test/*.cpp' ),
		target = 'test_bgen',
		use = 'bgen zlib',
		cxxflags = ['-std=c++11'],
		includes = 'test'
	)
