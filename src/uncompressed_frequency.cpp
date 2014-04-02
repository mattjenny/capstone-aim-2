#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>

using namespace std;

char HOMOZYGOUS = '1';
char HETEROZYGOUS = '2';

unsigned int get_size(string filename)
{
  ifstream f(filename);
  f.seekg(0, f.end);
  int size = f.tellg();
  return size;
}

unsigned int rank_query(string filename, char variant, int index)
{
  ifstream f(filename);
  if(!f)
    {
      cout << "Could not open file, enter valid file." << endl;
      return 0;
    }
  int count = 0;
  for(int i=0; i <= index; i++)
    {
      char v = f.get();
      if(v == variant)
	count++;
    }
  return count;
}

string get_high_frequency_variants(string filename, unsigned int population_size, unsigned int index, double freq_cutoff) {
	
	// Bounds check
	if (get_size(filename) < population_size*index) return "";

	stringstream ss;
	unsigned int cumulative_rank = 0;
	for (unsigned int i=0; i<index; i++) {
		unsigned int homozygous_count = rank_query(filename, HOMOZYGOUS, (i+1)*population_size - 1);
		unsigned int heterozygous_count = rank_query(filename, HETEROZYGOUS, (i+1)*population_size - 1);
		unsigned int row_count = homozygous_count + heterozygous_count - cumulative_rank;
		if (((double) row_count)/population_size >= freq_cutoff) 
		  {
		    ss << i;
		    ss << " ";
		  }
		cumulative_rank = homozygous_count + heterozygous_count;
	}

	return ss.str();
}

int main(int argc, char** argv) {
  if(argc < 2)
    {
      cout << "Usage: " << argv[0] << " filename frequency" << endl;
    }
 string filename = argv[1];
 double frequency = atof(argv[2]);
  
 string s = get_high_frequency_variants(filename, 20, 10, frequency);
 cout << "RESULT STRING: (row indices where percent variants >= cutoff): " << s << endl;
}
