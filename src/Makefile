CXX	= g++
CFLAGS = -std=c++0x
OBJS = BitPrinter.o BitReader.o NavarroSeq.o capstone_query.o
OBJS2 = BitPrinter.o BitReader.o NavarroSeq.o navarro_compress.o
PROG = CapstoneQuery 
PROG2 = NavarroCompress
.SUFFIXES: .o .cpp

all: $(PROG) $(PROG2)

.cpp.o:
	$(CXX) $(CFLAGS) -c -o $@ $<

$(PROG): $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS)

$(PROG2): $(OBJS2)
	$(CXX) $(CFLAGS) -o $@ $(OBJS2)

clean:
	/bin/rm -f *.o *~

