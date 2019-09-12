FLAGS = -g -std=c++11 -lz \
-I genfile/include \
-I db/include \
-I 3rd_party/boost_1_55_0 \
-I 3rd_party/zstd-1.1.0 \
-I 3rd_party/zstd-1.1.0/lib \
-I 3rd_party/zstd-1.1.0/lib/common \
-I 3rd_party/zstd-1.1.0/lib/compress \
-I 3rd_party/zstd-1.1.0/lib/decompress \
-I 3rd_party/sqlite3 \
-I include/3rd_party/sqlite3 \
-D SQLITE_ENABLE_COLUMN_METADATA \
-D SQLITE_ENABLE_STAT4 \
-D SQLITE_MAX_EXPR_DEPTH=10000 \
-D SQLITE_USE_URI=1 \
-Wno-unused-local-typedefs \
-fPIC -O3

bgen_to_vcf: example/bgen_to_vcf.cpp $(wildcard 3rd_party/zstd-1.1.0/*.o) $(wildcard 3rd_party/sqlite3/*.o)
	g++ ${FLAGS} -o build/bgen_to_vcf example/bgen_to_vcf.cpp src/*.cpp
