#include "NavarroSeq.h"
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

char HOMOZYGOUS = 0x31;
char HETEROZYGOUS = 0x32;

string get_high_frequency_variants(string filename, unsigned int population_size, unsigned int numrows, double freq_cutoff) {
	
	// Bounds check
	if (NavarroSeq::get_size(filename) < population_size*numrows) return "";

	stringstream ss;
	unsigned int cumulative_rank = 0;
	for (unsigned int i=0; i<numrows; i++) {
		unsigned int homozygous_count = NavarroSeq::rank(filename, HOMOZYGOUS, (i+1)*population_size - 1);
		unsigned int heterozygous_count = NavarroSeq::rank(filename, HETEROZYGOUS, (i+1)*population_size - 1);
		unsigned int row_count = homozygous_count + heterozygous_count - cumulative_rank;
		if (((double) row_count)/population_size >= freq_cutoff) ss << i;
		cumulative_rank = homozygous_count + heterozygous_count;
	}

	return ss.str();
}

int main(int argc, char** argv) {

	bool compress = false;
	bool decompress = false;
	string fileIn = "";
	string fileOut = "";

	if (argc > 1) {
		if (strcmp("-compress", argv[1]) == 0) {
			compress = true;
			fileIn = argv[2];
			fileOut = argv[3];
		} 
		else if (strcmp("-decompress", argv[1]) == 0) {
			decompress = true;
			fileIn = argv[2];
		} else {
			fileIn = argv[1];
		}
	} else {
		cout << "Not enough arguments!" << endl;
		return 0;
	}

	if (compress) {
		NavarroSeq::compress(fileIn, fileOut);
	} else if (decompress) {
		NavarroSeq::decompress(fileIn);
	} else {
		string s = get_high_frequency_variants(fileIn, 20, 10, 0.5);
		cout << "RESULT STRING: (row indices where percent variants >= cutoff): " << s << endl;
	}
}
