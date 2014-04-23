#include <fstream>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <time.h>   

using namespace std;

int main(int c, char** argv) {

	if (c < 3) return 0;
	unsigned long filesize = atol(argv[1]);
	string fname(argv[2]);

	std::ofstream ofs(fname, std::ofstream::out);

	srand(time(NULL));

	for (size_t i = 0; i<filesize; i++) {
		char c = (char) (rand() % 4 + 0x30);
		ofs.put(c);
	}

	ofs.close();

	return 0;
}