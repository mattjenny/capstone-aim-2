CXX	= g++
.SUFFIXES: .o .cpp

all: NavarroSeq

BitPrinter.o: BitPrinter.cpp BitPrinter.h
	$(CXX) -c BitPrinter.cpp

BitReader.o: BitReader.cpp BitReader.h
	$(CXX) -c BitReader.cpp

NavarroSeq: NavarroSeq.cpp BitReader.o BitPrinter.o ETable.h
	$(CXX) -o NavSeq BitReader.o BitPrinter.o NavarroSeq.cpp

clean:
	/bin/rm -f *.o *~

