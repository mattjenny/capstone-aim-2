#include "NavarroSeq.h"
#include <string>
#include <math.h>
#include <unordered_set>
#include <stdio.h>
#include <vector>
#include <list>
#include <limits>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

size_t nseq_n, nseq_r, nseq_u;
vector<char> nseq_alphabet;

union {
	unsigned int integer;
	char byte[4];
} four_byte_union;

void print_int(std::ofstream & ofs, unsigned int integer) {
	four_byte_union.integer = integer;
	ofs.put(four_byte_union.byte[0]);
	ofs.put(four_byte_union.byte[1]);
	ofs.put(four_byte_union.byte[2]);
	ofs.put(four_byte_union.byte[3]);
}

unsigned int get_int(std::ifstream & ifs) {
	char bytes[4];
	bytes[0] = ifs.get();
	bytes[1] = ifs.get();
	bytes[2] = ifs.get();
	bytes[3] = ifs.get();
	four_byte_union.byte[0] = bytes[0];
	four_byte_union.byte[1] = bytes[1];
	four_byte_union.byte[2] = bytes[2];
	four_byte_union.byte[3] = bytes[3];

	return four_byte_union.integer;
}

// Private constructor: instances used internally only
NavarroSeq::NavarroSeq(size_t n, size_t r, size_t u, list<char> alphabet) 
{
	nseq_n = n;
	nseq_r = r;
	nseq_u = u;
	for (list<char>::iterator it = alphabet.begin(); it != alphabet.end(); ++it) {
		nseq_alphabet.push_back(*it);
	}
}

void NavarroSeq::compress(string in_fname, string out_fname)
{
	/*
	* n is the length of the full sequence s
	* r is the number of characters in the alphabet
	* u is the block size
	* num_blocks is the number of u-length blocks in the sequence
	* i,j,k are counters used throughout this method
	*/
	size_t n,r,u,num_blocks,i,j,k;

	std::ifstream ifs (in_fname);

	// Get file size:
	ifs.seekg(0, ifs.end);
	n = ifs.tellg();
	ifs.seekg(0, ifs.beg);

	// Figure out the alphabet that this sequence uses
	unordered_set<char> unique_chars;
	for (i=0; i<n; i++) {
		unique_chars.insert(ifs.get());
	}
	ifs.seekg(0, ifs.beg);
	r = unique_chars.size();

	// Sort the alphabet
	list<char> alphabet;
	for(unordered_set<char>::iterator set_it = unique_chars.begin(); set_it != unique_chars.end(); ++set_it) {
		alphabet.push_back(*set_it);
	}
	alphabet.sort();

	// Calculate block length u and number of blocks
	u = floor(0.5*log((double) n)/log((double) r));
	if (u <= 0) u = 1; //block size must be positive
	num_blocks = floor((n+0.0)/u);

	// Create private instance of NavarroSeq
	NavarroSeq* nseq = new NavarroSeq(n, r, u, alphabet);
	
	// Enumerate and index all possible COMBINATIONS of alphabet characters
	vector<vector<unsigned int> > combinations = get_all_combinations(unique_chars.size(), u);

	//Create E tables:
	vector<E_table*> E_table_ptrs; // Mapping from R index to E table
	for (i=0; i<combinations.size(); i++) {
		// Create E table for R[i]
		E_table* etable = nseq->get_etable(combinations.at(i));
		E_table_ptrs.push_back(etable);
	}

	/*
	* r_vals: the fixed-length Ri identifiers for the combination of each block i
	* i_vals: the variable-length Ii identifiers for the permutation of each block i
	* l_partial_sums: the sum of I encoding lengths up to block i
	* n_partial_sums: the sum of character counts for each character up to block i
	*/
	vector<unsigned int> r_vals (num_blocks, 0); //Note: using unsigned int wastes space!
	vector<unsigned int> i_vals (num_blocks, 0); //Note: using unsigned int wastes space!
	vector<unsigned int> l_partial_sums (num_blocks+1, 0); //NOT exactly the definiteion from Navarro; num G entries
	vector<unsigned int> n_partial_sums ((num_blocks+1)*r, 0); // 2-D: blocks and chars

	l_partial_sums.push_back(0);

	unsigned int curr_r_index, curr_i_index, curr_l;
	for (i=0; i<num_blocks; i++) { //Iterate over every full block
		//current_block = s.substr(i*u, u);
		char data[u+1];
		ifs.read(data, u);
		data[u] = 0x00;
		string current_block (data);

		// Inefficient: we don't really care about compression efficiency, but this could be better
		// Figure out which index in R this is
		vector<unsigned int> curr_combination (r,0); // Tranform string into character vector
		for (j=0; j<current_block.length(); j++) {
			for (k=0; k<r; k++) {
				if (current_block.at(j) == nseq_alphabet.at(k)) {
					curr_combination.at(k)++;
					break;
				}
			}
		}
		for (j=0; j<combinations.size(); j++) {
			bool is_match = true;
			for (k=0; k<r; k++) {
				if (curr_combination.at(k) != combinations.at(j).at(k)) {
					is_match = false;
					break;
				}
			}
			if (is_match) {
				r_vals.at(i) = j;
				curr_r_index = j;
				break;
			}
		}

		// Next, find permutation Ii
		E_table* table = E_table_ptrs.at(curr_r_index);
		curr_l = ceil(log(table->entries.size()));
		l_partial_sums.at(i+1) = (l_partial_sums.at(i) + curr_l);

		for (j=0; j<table->entries.size(); j++) {
			G_entry* entry = table->entries.at(j);
			if(current_block.compare(entry->sequence) == 0) {
				i_vals.at(i) = j;
				// update n_partial_sums:
				for (k=0; k<r; k++) {
					n_partial_sums.at((i+1)*r+k) = n_partial_sums.at(i*r+k) + curr_combination.at(k);
				}
				break;
			}
		}
	}

	// Now for the remainder that didn't fit into a block:
	int rmdr_length = n%u;
	char data[rmdr_length+1];
	ifs.read(data, rmdr_length);
	data[rmdr_length] = 0x00;
	string rmdr (data);

	std::ofstream ofs(out_fname);

	print_int(ofs, n);
	print_int(ofs, r);
	print_int(ofs, u);
	print_int(ofs, num_blocks);
	print_int(ofs, combinations.size());

	for (list<char>::iterator alphabet_it = alphabet.begin(); alphabet_it != alphabet.end(); ++alphabet_it) {
		ofs.put(*alphabet_it);
	}

	for (i=0; i<r_vals.size(); i++) {
		print_int(ofs, r_vals.at(i));
	}
	
	for (i=0; i<i_vals.size(); i++) {
	  	print_int(ofs, i_vals.at(i));
	}

	for (i=0; i<l_partial_sums.size(); i++) {
	  	print_int(ofs, l_partial_sums.at(i));
	}

	for (i=0; i<n_partial_sums.size(); i++) {
		print_int(ofs, n_partial_sums.at(i));
	}

	unsigned int total_g_table_depth = 0;
	print_int (ofs, total_g_table_depth);
	for (i=0; i<E_table_ptrs.size(); i++) {
		E_table* etable = E_table_ptrs.at(i);
		total_g_table_depth += etable->entries.size();
		print_int(ofs, total_g_table_depth);
	}

	for (i=0; i<E_table_ptrs.size(); i++) {
		E_table* etable = E_table_ptrs.at(i);
		for (unsigned int j=0; j<etable->entries.size(); j++) {
			G_entry* entry = etable->entries.at(j);
			// Output compressed code
			ofs << entry->sequence;
			for(unsigned int k=0; k<entry->ranks.size(); k++) {
				  ofs.put(entry->ranks.at(k));
			}
		}
	}

	ofs << rmdr;
}

string NavarroSeq::decompress(string filename)
{

	std::ifstream ifs (filename, std::ifstream::in);

	unsigned int i,j,k;
	unsigned int n = get_int(ifs);
	unsigned int r = get_int(ifs);
	unsigned int u = get_int(ifs);
	unsigned int num_blocks = get_int(ifs);
	unsigned int num_combos = get_int(ifs);

	// Skip stored alphabet
	for (i=0; i<r; i++) {
		ifs.get();
	}

	vector<unsigned int> r_vals (num_blocks, 0); //Note: using unsigned int wastes space!
	vector<unsigned int> i_vals (num_blocks, 0); //Note: using unsigned int wastes space!
	vector<unsigned int> l_partial_sums (num_blocks+1, 0);
	vector<unsigned int> E_table_depths (num_combos+1);
	vector<E_table*> E_table_ptrs; // Mapping from R index to E table

	// Populate r_vals
	for (i=0; i<num_blocks; i++) {
		r_vals.at(i) = get_int(ifs);
	}

	// Populate i_vals
	for (i=0; i<num_blocks; i++) {
		i_vals.at(i) = get_int(ifs);
	}

	// Populate l_partial_sums
	for (i=0; i<num_blocks+1; i++) {
		l_partial_sums.at(i) = get_int(ifs);
	}

	// Skip n_partial_sums
	ifs.ignore((num_blocks+1)*r*4+4, EOF); 

	// Populate E Table Depths
	for (i=0; i<num_combos+1; i++) {
		E_table_depths.at(i) = get_int(ifs);
	}

	// Populate E Tables
	for (i=0; i<num_combos; i++) {
		E_table* table = new E_table();
		unsigned int num_entries = E_table_depths.at(i+1) - E_table_depths.at(i);
		for (j=0; j<num_entries; j++) {
			char c_sequence[u+1];
			for (k=0; k<u; k++) {
				c_sequence[k] = ifs.get();
			}
			c_sequence[u] = 0x00;
			string sequence (c_sequence);
			vector<char> ranks (r*u, 0);
			for (k=0; k<r*u; k++) {
				ranks.at(k) = ifs.get();
			}
			G_entry* entry = new G_entry(sequence, ranks);
			table->add_entry(entry);
		}
		E_table_ptrs.push_back(table);
	}

	char c;
	string rmdr = "";
	while ((c = ifs.get()) != EOF) {
		rmdr = rmdr + c;
	}

	string s = "";
	
	for (unsigned int i=0; i<num_blocks; i++) {
		E_table* table = E_table_ptrs.at(r_vals.at(i));
		G_entry* entry = table->entries.at(i_vals.at(i));
		s += entry->sequence;
	}

	cout << "DECOMPRESSED STRING = " << s << rmdr << endl;
	
	return s;
}

char NavarroSeq::access(string fname, int index)
{
	std::ifstream ifs (fname, std::ifstream::in);

	unsigned int i,j,k;
	unsigned int n = get_int(ifs);
	unsigned int r = get_int(ifs);
	unsigned int u = get_int(ifs);
	unsigned int num_blocks = get_int(ifs);
	unsigned int num_combos = get_int(ifs);

	// Skip stored alphabet
	for (i=0; i<r; i++) {
		ifs.get();
	}

	unsigned int block = floor(index/u);
	unsigned int l = index - block*u;

	ifs.seekg((5 + block) * 4 + r);
	unsigned int r_val = get_int(ifs);

	ifs.seekg((5 + num_blocks + block) * 4 + r);
	unsigned int i_val = get_int(ifs);

	ifs.seekg((7+num_blocks*(3+r) + r + r_val)*4 + r);
	unsigned int E_depth = get_int(ifs);

	ifs.seekg((8+num_blocks*(3+r)+r+num_combos)*4 + E_depth*((r + 1) * u) + l + r);
	char retval = ifs.get();

	return retval;
}

unsigned int NavarroSeq::rank(string fname, char c, int index)
{

	std::ifstream ifs (fname, std::ifstream::in);

	unsigned int char_index_in_alphabet = 0; // hard coded for now

	unsigned int i,j,k;
	unsigned int n = get_int(ifs);
	unsigned int r = get_int(ifs);
	unsigned int u = get_int(ifs);
	unsigned int num_blocks = get_int(ifs);
	unsigned int num_combos = get_int(ifs);

	// Find in stored alphabet
	for (i=0; i<r; i++) {
		char current_char = ifs.get();
		if (current_char == c) char_index_in_alphabet = i;
	}

	unsigned int block = floor(index/u);
	unsigned int l = index - block*u;

	ifs.seekg((5 + block) * 4 + r);
	unsigned int r_val = get_int(ifs);

	ifs.seekg((5 + num_blocks + block) * 4 + r);
	unsigned int i_val = get_int(ifs);

	ifs.seekg((7 + 3*num_blocks + r*block + char_index_in_alphabet) * 4 + r);
	unsigned int n_partial_sum = get_int(ifs);

	ifs.seekg((7+num_blocks*(3+r)+r + r_val)*4 + r);
	unsigned int E_depth = get_int(ifs);

	ifs.seekg((8+num_blocks*(3+r)+r+num_combos)*4 + E_depth*((r + 1) * u) + u + u*char_index_in_alphabet + l + r);
	char block_rank = ifs.get();

	unsigned int rank = n_partial_sum + block_rank;

	//cout << "RANK for index " << index << " = " << n_partial_sum << " + " << (int)block_rank << " = " << rank << endl;

	return rank;

}

unsigned int NavarroSeq::select(string fname, char c, int index)
{
	return 55;
}

unsigned int NavarroSeq::get_size(string fname) {
	std::ifstream ifs (fname, std::ifstream::in);
	return get_int(ifs);
}

vector<vector<unsigned int> > NavarroSeq::get_all_combinations(size_t r, size_t u) // r is alphabet size; u is block size
{
	unsigned int i;
	vector<vector<unsigned int> > combinations; // store a list of all combinations

	for (i = 0; i<r; i++)
	{
		/*
		* Create a vector of counts for each character
		* Vector of length r
		* values[0] has counts for alphabet[0], values[1] for alphabet[1]...
		*/
		vector<unsigned int> values (r,0);
		values.at(i)++; // Set initial count to current character
		enumerate(values, combinations, i, r, u-1); // add all enumerations with this prefix
	}
	return combinations;
}

void NavarroSeq::enumerate(vector<unsigned int>& values, vector<vector<unsigned int> >& combos, unsigned int min, unsigned int r, unsigned int u)
{
	if (u <= 0) //We already have a sequence of length u
	{
		combos.push_back(values); // Add it to the list
		return;
	}
	else
	{
		unsigned int i;
		for (i = min; i<r; i++) // iterate through all possible characters >= prefix
		{
			vector<unsigned int> values_copy (values);
			values_copy.at(i)++;
			enumerate(values_copy, combos, i, r, u-1);
		}
	}
}

void NavarroSeq::get_etable_rows(string prefix, vector<char> ranks, vector<unsigned int> combination, E_table* table)
{
	unsigned int i,j,combination_sum = 0;
	for (i=0; i<combination.size(); i++) {
		combination_sum += combination.at(i);
	}
	if (combination_sum == 0) {
		G_entry* entry = new G_entry(prefix, ranks);
		table->add_entry(entry);
	}
	
	for (i=0; i<combination.size(); i++) {
		if(combination.at(i) == 0) continue;
		string next_prefix = prefix + nseq_alphabet.at(i);
		vector<unsigned int> remaining_combo (combination);
		remaining_combo.at(i)--;
		for (j=i*nseq_u + prefix.length(); j<(i+1)*nseq_u; j++) {
			ranks.at(j)++;
		}
		get_etable_rows(next_prefix, ranks, remaining_combo, table);
	}
}

E_table* NavarroSeq::get_etable(vector<unsigned int> combination)
{
	E_table* table = new E_table();
	unsigned int i,j;
	for (i=0; i<combination.size(); i++) {
		if(combination.at(i) == 0) continue;
		string prefix = "";
		prefix += nseq_alphabet.at(i);
		vector<char> ranks (nseq_r*nseq_u, 0);
		vector<unsigned int> remaining_combo (combination);
		remaining_combo.at(i)--;
		for (j=i*nseq_u; j<(i+1)*nseq_u; j++) {
			ranks.at(j)++;
		}
		get_etable_rows(prefix, ranks, remaining_combo, table);
	}
	return table;
}

/*int main() {
	// test sequence
	//NavarroSeq::compress("testin.txt", "testout.txt");
	//NavarroSeq::compress("abbabaababbababcabbabaababbababcabbabaababbababcabbabaababbababcabbabaababbababcabc");
	//NavarroSeq::decompress("testout.txt");
	//cout << NavarroSeq::access("testout.txt",31) << endl;
	NavarroSeq::rank("testout.txt",0x62,31);

	return 0;
}*/
