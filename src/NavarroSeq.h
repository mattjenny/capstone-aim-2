# ifndef NAVARROSEQ_H
# define NAVARROSEQ_H

#include <string>
#include <vector>
#include <unordered_map>

using namespace std;

class NavarroSeq
{
public:
	static string compress(string s);
	static string decompress(string s);
	static char access(string s, int index);
	static unsigned int rank(string s, char c, int index);
	static unsigned int select(string s, char c, int index);

private:
	static vector<vector<unsigned int> > get_all_combinations(size_t r, size_t u);
	static void enumerate(vector<unsigned int>& values, vector<vector<unsigned int> >& combos, unsigned int min, unsigned int r, unsigned int u);
};

# endif