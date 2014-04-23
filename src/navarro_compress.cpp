#include "NavarroSeq.h"

using namespace std;

int main(int argc, char** argv) {
	if (argc < 3) return 0;
	NavarroSeq::compress(argv[1],  argv[2]);
}