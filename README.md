#Matt Jenny and Linda Zhang

##Repository for Implementation code for Aim 2

###Goal 1: Implement Navarro's multiary wavelet trees
Compress, decompress, access, and rank are implemented as public static methods in NavarroSeq.cpp.  frequency_query.cpp contains code to solve the following query: "given a set of individuals and a set of possible variants, find all variants such that at least some given percentage of the individuals have the mutation."  The main method in frequency_compress.cpp currently compresses the contents of "testin.txt" to "testout.txt" with the -compress flag, and otherwise runs the aforementioned query on testout.txt.

###Goal 2: Implement Barbay structure