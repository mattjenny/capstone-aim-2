#include "NavarroSeq.h"
#include <string.h>
#include <fstream>
#include <sstream>

using namespace std;

char HOMOZYGOUS = 0x31;
char HETEROZYGOUS = 0x32;

string get_persons_with_variant(string filename, char variant, unsigned int index, unsigned int population_size)
{
  stringstream persons;
  unsigned int start_location = population_size*index;
  unsigned int end_location = start_location+population_size;

  for(unsigned int i=start_location; i<=end_location; i++)
    {
      char v = NavarroSeq::access(filename, i);
      if(v == variant)
	{
	  persons << i-start_location;
	  persons << " ";
	}
    }
  return persons.str();
}

int main(int argc, char** argv)
{
  if(argc < 4)
    {
      cout << "Usage: " << argv[0] << " filename variant index populationSize" << endl;
    }
  string filename = argv[1];
  char variant;
  unsigned int index = atoi(argv[3]);
  unsigned int populationSize = atoi(argv[4]);

  if(strcmp(argv[2], "1") == 0)
    {
      variant = HOMOZYGOUS;
      cout << "Homozygous" << endl;
    }
  else if(strcmp(argv[2], "2") == 0)
    {
      variant = HETEROZYGOUS;
      cout << "Heterozygous" << endl;
    }
  string s = get_persons_with_variant(filename, variant, index, populationSize);
  cout << "RESULT: " << s << endl;
}
