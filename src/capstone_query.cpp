#include "NavarroSeq.h"
#include <string.h>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

unsigned char HOMOZYGOUS = 0x31;
unsigned char HETEROZYGOUS = 0x32;

bool is_variant(unsigned char c) {
  return c == HOMOZYGOUS || c == HETEROZYGOUS;
}

string ca_query_1(string filename, unsigned int variant_index, unsigned int population_size)
{
  stringstream persons;
  unsigned int start_location = population_size*variant_index;
  unsigned int end_location = start_location+population_size;

  for(unsigned int i=start_location; i<end_location; i++)
    {
      unsigned char v = NavarroSeq::access(filename, i);
      if(is_variant(v)) 
      {
    	  persons << i-start_location;
    	  persons << " ";
    	}
    }
  return persons.str();
}

string uncompressed_query_1(string filename, unsigned int variant_index, unsigned int population_size) 
{
  stringstream persons;
  unsigned int start_location = population_size*variant_index;
  unsigned int end_location = start_location+population_size;
  std::ifstream ifs(filename);
  ifs.seekg(population_size*variant_index);
  for(unsigned int i=start_location; i<end_location; i++)
    {
      unsigned char v = (unsigned char) ifs.get();
      if(is_variant(v)) 
      {
        persons << i-start_location;
        persons << " ";
      }
    }
    ifs.close();
  return persons.str();
}

string ca_query_2(string filename, unsigned int num_variants, unsigned int population_size, double freq_cutoff) 
{
  // Bounds check
  if (NavarroSeq::get_size(filename) < population_size*num_variants) return "";

  stringstream ss;
  unsigned int cumulative_rank = 0;
  for (unsigned int i=0; i<num_variants; i++) {
    unsigned int homozygous_count = NavarroSeq::rank(filename, (i+1)*population_size - 1, HOMOZYGOUS);
    unsigned int heterozygous_count = NavarroSeq::rank(filename, (i+1)*population_size - 1, HETEROZYGOUS);
    unsigned int row_count = homozygous_count + heterozygous_count - cumulative_rank;
    if (((double) row_count)/population_size >= freq_cutoff) {
      ss << i << " ";
    }
    cumulative_rank = homozygous_count + heterozygous_count;
  }

  return ss.str();
}

string uncompressed_query_2(string filename, unsigned int num_variants, unsigned int population_size, double freq_cutoff) 
{
  stringstream ss;
  std::ifstream ifs(filename);
  unsigned int row_count;
  unsigned char c;
  for (unsigned int i=0; i<num_variants; i++) {
    row_count = 0;
    for (unsigned int j=0; j<population_size; j++) {
      c = (unsigned char) ifs.get();
      if (is_variant(c)) row_count++;
    }
    if (((double) row_count)/population_size >= freq_cutoff) {
      ss << i << " ";
    }
  }
  ifs.close();
  return ss.str();
}

string run_query(string query_arg, string filename, unsigned int variant_arg, unsigned int population_size)
{
  if (query_arg.compare("-c1") == 0) {
    return ca_query_1(filename, variant_arg, population_size);
  } 
  else if (query_arg.compare("-u1") == 0) {
    return uncompressed_query_1(filename, variant_arg, population_size);
  }
  else {
    return "Invalid flag. Type '-help' for more.";
  }
}

int main(int argc, char** argv)
{
  if (argc > 1) {
    if (strcmp(argv[1], "-help") == 0) {
      cout << "Capstone Query V1.0" << endl;
      cout << "@author Matt Jenny (mvj5fs)" << endl << endl;
      cout << "Usage: " << endl;
      cout << "       For query 1: " << argv[0] << " $(FLAG) [your_file_name] [variant_index] [population_size]" << endl;
      cout << "       For query 2: " << argv[0] << " $(FLAG) [your_file_name] [variant_count] [population_size] [frequency cutoff]" << endl << endl;
      cout << "Flags: " << endl;
      cout << "       -c1 : Compression-aware algorithm for query 1" << endl;
      cout << "       -u1 : Uncompressed algorithm for query 1" << endl;
      cout << "       -c2 : Compression-aware algorithm for query 2" << endl;
      cout << "       -u2 : Uncompressed algorithm for query 2" << endl;
      return 0;
    }
  } 
  if (argc < 5)
  {
      cout << "Wrong number of arguments; type '-help' for more." << endl;
      return 0;
  }
  string query_arg = string(argv[1]);
  string filename = argv[2];
  unsigned int variant_arg = atoi(argv[3]);
  unsigned int population_size = atoi(argv[4]);
  string s;
  if (query_arg.compare("-c1") == 0) {
    s = ca_query_1(filename, variant_arg, population_size);
  } 
  else if (query_arg.compare("-u1") == 0) {
    s = uncompressed_query_1(filename, variant_arg, population_size);
  }
  else if (argc > 5) {
    double freq_cutoff = atof(argv[5]);
    if (query_arg.compare("-c2") == 0) {
      s = ca_query_2(filename, variant_arg, population_size, freq_cutoff);
    } 
    else if (query_arg.compare("-u2") == 0) {
      s = uncompressed_query_2(filename, variant_arg, population_size, freq_cutoff);
    } else {
      s = "Invalid flag. Type '-help' for more.";
    }
  } else {
    s = "Invalid flag or too few arguments. Type '-help' for more.";
  }
  cout << s << endl;
  return 0;
}
