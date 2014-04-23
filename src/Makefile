CXX	= g++
CFLAGS = -std=c++0x
OBJS = BitPrinter.o BitReader.o NavarroSeq.o capstone_query.o
PROG = CapstoneQuery
.SUFFIXES: .o .cpp

all: $(PROG)

.cpp.o:
	$(CXX) $(CFLAGS) -c -o $@ $<

$(PROG): $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS)

clean:
	/bin/rm -f *.o *~

