# ifndef NAVARROSEQ_H
# define NAVARROSEQ_H

#include <string>
#include <vector>
#include <unordered_map>
#include "ETable.h"
#include <list>
#include <limits>
#include <iostream>

using namespace std;

class NavarroSeq
{
public:
	static void compress(string in_fname, string out_fname);
	static string decompress(string filename);
	static char access(string fname, int index);
	static unsigned int rank(string fname, char c, int index);
	static unsigned int select(string fname, char c, int index);

private:
	NavarroSeq(size_t n, size_t r, size_t u, list<char> alphabet);
	static vector<vector<unsigned int> > get_all_combinations(size_t r, size_t u);
	static void enumerate(vector<unsigned int>& values, vector<vector<unsigned int> >& combos, unsigned int min, unsigned int r, unsigned int u);
	void get_etable_rows(string prefix, vector<char> ranks, vector<unsigned int> combination, E_table* table);
	E_table* get_etable(vector<unsigned int> combination);
};

# endif
