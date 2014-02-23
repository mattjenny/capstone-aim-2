# ifndef NAVARROSEQ_H
# define NAVARROSEQ_H

#include <string>

using namespace std;

class NavarroSeq
{
	static string compress(string s);
	static string decompress(string s);
	static char access(string s, int index);
	static unsigned int rank(string s, char c, int index);
	static unsigned int select(string s, char c, int index);
};

# endif