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

int main() {
	//NavarroSeq::compress("testin.txt", "testout.txt");
	string s = get_high_frequency_variants("testout.txt", 4, 4, 0.5);
	cout << "RESULT STRING: (row indices where percent variants > cutoff): " << s << endl;
}