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
	static string toBinStr(unsigned int val);
	static string compress(string s);
	static string decompress(vector<E_table*>& E_table_ptrs, vector<unsigned int>& r_vals, vector<unsigned int>& i_vals, vector<unsigned int>& l_partial_sums, vector<unsigned int>& n_partial_sums, string rmdr);
	static char access(unsigned int u, vector<E_table*>& E_table_ptrs, vector<unsigned int>& r_vals, vector<unsigned int>& i_vals, vector<unsigned int>& l_partial_sums, vector<unsigned int>& n_partial_sums, string rmdr, int index);
	static unsigned int rank(string s, char c, int index);
	static unsigned int select(string s, char c, int index);

private:
	NavarroSeq(size_t n, size_t r, size_t u, list<char> alphabet);
	static vector<vector<unsigned int> > get_all_combinations(size_t r, size_t u);
	static void enumerate(vector<unsigned int>& values, vector<vector<unsigned int> >& combos, unsigned int min, unsigned int r, unsigned int u);
	void get_etable_rows(string prefix, vector<char> ranks, vector<unsigned int> combination, E_table* table);
	E_table* get_etable(vector<unsigned int> combination);
};

# endif
