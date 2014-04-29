#include "NavarroSeq.h"
#include "BitPrinter.h"
#include "BitReader.h"
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

int BYTE_SIZE = 8;
int NUM_BYTES_IN_INT = 4;

/*
* n is the length of the full sequence s
* r is the number of characters in the alphabet
* u is the block size
* num_blocks is the number of u-length blocks in the sequence
*/
size_t n,r,u,num_blocks;
vector<char> alphabet;
vector<vector<unsigned int> > combinations;
vector<E_table*> E_table_ptrs; // Mapping from R index to E table
vector<unsigned int> curr_combination(r,0);
unsigned int large_block_size;


//Element bit sizes
unsigned int r_int_size;
unsigned int full_l_int_size;
unsigned int relative_l_int_size;
unsigned int n_int_size;

union {
	unsigned int integer;
	char byte[4];
} four_byte_union;

union {
	unsigned long thelong;
	char byte[8];
} eight_byte_union;

typedef struct {
	unsigned int value;
	unsigned int bits;
} L_val;

unsigned char on_deck = 0x00;
int num_bits_on_deck = 0;

/**********************************************************************************************************************************
**************************************** HELPER METHODS TO HANDLE BITWISE I/O *****************************************************
**********************************************************************************************************************************/


void print_char(std::ofstream & ofs, char c_in, unsigned int num_bits) {
	unsigned char c = (unsigned char) c_in;
	if (num_bits == 0) return;
	if (num_bits_on_deck + num_bits < BYTE_SIZE) {
		unsigned char added_bits = c << (8-num_bits-num_bits_on_deck);
		on_deck += added_bits;
		num_bits_on_deck += num_bits;
	}
	else {
		int shift_length = num_bits_on_deck + num_bits - 8;
		unsigned char shifted_c_filter = (unsigned char)0xff >> (8-(num_bits - shift_length));
		unsigned char shifted_c = (c >> shift_length) & shifted_c_filter;
		unsigned char retval = on_deck + (shifted_c);
		num_bits_on_deck = shift_length;
		on_deck = c << (8 - shift_length);
		ofs.put(retval);
	}
}

void print_int(std::ofstream & ofs, unsigned int integer, unsigned int num_bits) {

	if (num_bits > 32) {
		cout << "There aren't that many bits in an int." << endl;
		return;
	} else if (num_bits == 0) {
		return;
	}

	four_byte_union.integer = integer;

	char prev = on_deck;
	char next = 0;
	int i;

	int prev_bits_left = num_bits_on_deck;
	int mod_bits = num_bits % 8;
	if (mod_bits == 0) mod_bits = 8;
	int max_byte = (num_bits - 1)/8;
	print_char(ofs, four_byte_union.byte[max_byte], mod_bits);

	for (i=max_byte - 1; i>=0; i--) {
		print_char(ofs, four_byte_union.byte[i], 8);
	} 
}

unsigned char get_char(std::ifstream & ifs, int num_bits) {
	if (num_bits > 8) {
		cout << "Char length must not exceed eight bits!" << endl;
		return 0;
	}

	unsigned char retval;

	if (num_bits <= num_bits_on_deck) {
		retval = on_deck >> (8 - num_bits);
		on_deck = on_deck << num_bits;
		num_bits_on_deck -= num_bits;
	}
	else {
		unsigned char c = (unsigned char) ifs.get();
		retval = (on_deck >> (8 - num_bits)) + (c >> (8 - num_bits + num_bits_on_deck));
		on_deck = c << (num_bits - num_bits_on_deck);
		num_bits_on_deck = 8 - num_bits + num_bits_on_deck;
	}

	return retval;
}

unsigned int get_int(std::ifstream & ifs, int num_bits) {
	if (num_bits > 32) {
		cout << "There aren't " << num_bits << " bits in an int." << endl;
		return 0;
	} else if (num_bits == 0) {
		return 0;
	}

	for (int i=0; i<NUM_BYTES_IN_INT; i++) {
		if ((3-i)*8 >= num_bits) {
			four_byte_union.byte[3-i] = 0x00;
		}
		else if ((4-i)*8 <= num_bits) {
			four_byte_union.byte[3-i] = get_char(ifs, BYTE_SIZE);
		} 
		else {
			four_byte_union.byte[3-i] = get_char(ifs, num_bits % BYTE_SIZE);
		}
	}

	return four_byte_union.integer;
}

void print_remainder(std::ofstream & ofs) {
	ofs.put(on_deck);
}

void flush_remainder() {
	on_deck = 0x00;
	num_bits_on_deck = 0;
}

/**********************************************************************************************************************************
**************************************** END HELPER METHODS TO HANDLE BITWISE I/O *************************************************
**********************************************************************************************************************************/

/**********************************************************************************************************************************
**************************************** Helper Methods Created Before 4.13.2014 **************************************************
**********************************************************************************************************************************/

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

void get_etable_rows(string prefix, vector<char> ranks, vector<unsigned int> combination, E_table* table)
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
		string next_prefix = prefix + alphabet.at(i);
		vector<unsigned int> remaining_combo (combination);
		vector<char> current_ranks (ranks);
		remaining_combo.at(i)--;
		for (j=i*u + prefix.length(); j<(i+1)*u; j++) {
			current_ranks.at(j)++;
		}
		get_etable_rows(next_prefix, current_ranks, remaining_combo, table);
	}
}

E_table* get_etable(vector<unsigned int> combination)
{
	E_table* table = new E_table();
	unsigned int i,j;
	for (i=0; i<combination.size(); i++) {
		if(combination.at(i) == 0) continue;
		string prefix = "";
		prefix += alphabet.at(i);
		vector<char> ranks (r*u, 0);
		vector<unsigned int> remaining_combo (combination);
		remaining_combo.at(i)--;
		for (j=i*u; j<(i+1)*u; j++) {
			ranks.at(j)++;
		}
		get_etable_rows(prefix, ranks, remaining_combo, table);
	}
	return table;
}

/**********************************************************************************************************************************
************************************ End Helper Methods Created Before 4.13.2014 **************************************************
**********************************************************************************************************************************/

/**********************************************************************************************************************************
**************************************** Compression Helper Methods 4.13.2014 *****************************************************
**********************************************************************************************************************************/

size_t get_filesize(std::ifstream & ifs) {
	ifs.seekg(0, ifs.end);
	size_t fsize = ifs.tellg();
	ifs.seekg(0, ifs.beg);
	return fsize;
}

unordered_set<char> get_unique_chars(std::ifstream & ifs) {
	unsigned int percent_complete = 0;
	unordered_set<char> unique_chars;
	for (size_t i=0; i<n; i++) {
		if ((((float)i)/n)*100 > percent_complete + 1) {
			percent_complete++;
		}
		unique_chars.insert(ifs.get());
	}
	ifs.seekg(0, ifs.beg);
	return unique_chars;
}

vector<char> get_sorted_alphabet(std::ifstream & ifs) {
	int percent_complete = 0;

	// Figure out the alphabet that this sequence uses
	unordered_set<char> unique_chars = get_unique_chars(ifs);

	// Copy unique characters to a list and sort
	list<char> alphabet;
	for(unordered_set<char>::iterator set_it = unique_chars.begin(); set_it != unique_chars.end(); ++set_it) {
		alphabet.push_back(*set_it);
	}
	alphabet.sort();

	// Now copy sorted list into a vector to allow indexed access.  Not efficient, but it gets the job done.
	vector<char> nseq_alphabet;
	for (list<char>::iterator it = alphabet.begin(); it != alphabet.end(); ++it) {
		nseq_alphabet.push_back(*it);
	}
	return nseq_alphabet;
}

string get_next_block(std::ifstream & ifs) {
	char data[u+1];
	ifs.read(data, u);
	data[u] = 0x00;
	string current_block (data);
	return current_block;
}

unsigned int get_r_index(string block) {
	// Figure out which index in R this is
	curr_combination.clear();
	curr_combination.resize(r,0);
	size_t j,k;
	for (j=0; j<block.length(); j++) {
		for (k=0; k<r; k++) {
			if (block.at(j) == alphabet.at(k)) {
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
			return j;
		}
	}
	return 0; // should add an error check here
}

unsigned int get_i_index(E_table* table, string block) {
	for (size_t j=0; j<table->entries.size(); j++) {
		G_entry* entry = table->entries.at(j);
		if(block.compare(entry->sequence) == 0) {
			return j;
		}
	}

	return 0; //should do error checking here
}

void NavarroSeq::print_header(std::ofstream & ofs) {
	print_int(ofs, n, 32);
	print_int(ofs, r, 32);
	print_int(ofs, u, 32);
	print_int(ofs, num_blocks, 32);
	print_int(ofs, combinations.size(), 32);

	for (vector<char>::iterator alphabet_it = alphabet.begin(); alphabet_it != alphabet.end(); ++alphabet_it) {
		print_char(ofs, *alphabet_it, 8);
	}
}

void NavarroSeq::print_E_table_info(std::ofstream & ofs) {
	unsigned int total_g_table_depth = 0;
	print_int (ofs, total_g_table_depth, 32);
	for (size_t i=0; i<E_table_ptrs.size(); i++) {
		E_table* etable = E_table_ptrs.at(i);
		total_g_table_depth += etable->entries.size();
		print_int(ofs, total_g_table_depth, 32);
	}

	unsigned int e_table_rank_size = ceil(log2(u+1));
	for (size_t i=0; i<E_table_ptrs.size(); i++) {
		E_table* etable = E_table_ptrs.at(i);
		for (unsigned int j=0; j<etable->entries.size(); j++) {
			G_entry* entry = etable->entries.at(j);
			for (unsigned int k=0; k<u; k++) {
				print_char(ofs, entry->sequence.at(k), 8);
			}
			for(unsigned int k=0; k<entry->ranks.size(); k++) {
				print_char(ofs, entry->ranks.at(k), e_table_rank_size);
			}
		}
	}
}

string NavarroSeq::parse_input_data(std::ifstream & ifs, string r_fname, string i_fname, string l_fname, string n_fname) {
	
	unsigned int percent_complete=0, current_r=0, current_i=0, current_l=0, current_l_size=0, current_i_size=0, prev_l=0, last_full_l_sum = 0, current_n=0;
	vector<unsigned int> prev_n(r, 0);
	vector<unsigned int> last_full_n_sum(r, 0);

	std::ofstream ofs_r_vals(r_fname, std::ofstream::out);
	std::ofstream ofs_i_vals(i_fname, std::ofstream::out);
	std::ofstream ofs_l_vals(l_fname, std::ofstream::out);
	std::ofstream ofs_n_vals(n_fname, std::ofstream::out);

	BitPrinter r_printer(&ofs_r_vals);
	BitPrinter i_printer(&ofs_i_vals);
	BitPrinter l_printer(&ofs_l_vals);
	BitPrinter n_printer(&ofs_n_vals);

	unsigned int r_bits = 0;
	unsigned int l_super_bits = 0;
	unsigned int l_rel_bits = 0;
	unsigned int i_bits = 0;
	unsigned int n_super_bits = 0;
	unsigned int n_rel_bits = 0;
	unsigned int num_superblocks = 0;
	unsigned int num_rel = 0;
	unsigned long last_n = 0;

	l_printer.print_int(0, full_l_int_size);
	for (size_t i=0; i<r; i++) {
		n_printer.print_int(0, full_l_int_size);
	}

	for (size_t i=0; i<num_blocks; i++) { //Iterate over every full block
		if ((((float)i)/num_blocks)*100 > percent_complete + 1) {
			percent_complete++;
		}

		string current_block = get_next_block(ifs);
		current_r = get_r_index(current_block);
		E_table* table = E_table_ptrs.at(current_r);
		current_i = get_i_index(table, current_block);
		current_i_size = ceil(log2(table->entries.size()));
		
		r_printer.print_int(current_r, r_int_size);
		i_printer.print_int(current_i, current_i_size);
		r_bits += r_int_size;
		i_bits += current_i_size;

		/*
		* L and N partial sums are a bit trickier.  If this sequence begins a new L and N 'superblock', we store the full partial sums.
		* Otherwise, we store just the relative sum since the beginning of the current 'superblock'.
		*/
		if ((i+1) % large_block_size == 0) {
			num_superblocks++;
			current_l = (prev_l + current_i_size + last_full_l_sum);
			l_printer.print_int(current_l, full_l_int_size);
			l_super_bits += full_l_int_size;
			for (size_t j=0; j<r; j++) {
				current_n = last_full_n_sum.at(j) + prev_n.at(j) + curr_combination.at(j);
				n_printer.print_int(current_n, full_l_int_size);
				n_super_bits += full_l_int_size;
				last_full_n_sum.at(j) = current_n;
				prev_n.at(j) = 0;
			}

			last_full_l_sum = current_l;
			prev_l = 0;
		} else {
			num_rel++;
			current_l = (prev_l + current_i_size);
			l_printer.print_int(current_l, relative_l_int_size);
			l_rel_bits += relative_l_int_size;
			for (size_t j=0; j<r; j++) {
				current_n = prev_n.at(j) + curr_combination.at(j);
				n_printer.print_int(current_n, relative_l_int_size);
				last_n = last_n << relative_l_int_size;
				last_n += current_n;
				n_rel_bits += relative_l_int_size;
				prev_n.at(j) = current_n;
			}
			
			prev_l = current_l;
		}
	}

	unsigned int rmdr_length = n%u;
	char data[rmdr_length+1];
	ifs.read(data, rmdr_length);
	data[rmdr_length] = 0x00;
	string rmdr (data);

	r_printer.print_remainder();
	i_printer.print_remainder();
	l_printer.print_remainder();
	n_printer.print_remainder();

	ofs_r_vals.close();
	ofs_i_vals.close();
	ofs_l_vals.close();
	ofs_n_vals.close();

	return rmdr;
}

void init(std::ifstream & ifs) {

	n = get_filesize(ifs);
	alphabet = get_sorted_alphabet(ifs);
	r = alphabet.size();
	u = floor(0.5*log2((double) n)/log2((double) r));
	
	if (u <= 0) u = 1; //block size must be positive
	num_blocks = floor(((double)n)/u);
	
	// Enumerate and index all possible COMBINATIONS of alphabet characters
	combinations = NavarroSeq::get_all_combinations(r, u);

	//Create E tables:
	for (size_t i=0; i<combinations.size(); i++) {
		// Create E table for R[i]
		E_table* etable = get_etable(combinations.at(i));
		E_table_ptrs.push_back(etable);
	}

	//Calculate bit sizes of each compressed element:
	r_int_size = ceil(log2(combinations.size()));
	large_block_size = ceil(log2(n*log2(r)));
	full_l_int_size = ceil(log2(n*log2(r)));
	relative_l_int_size = ceil(log2(u * ceil(log2(r)) * ceil(log2(n*log2(r)))));
}

/**********************************************************************************************************************************
************************************ End Compression Helper Methods 4.13.2014 *****************************************************
**********************************************************************************************************************************/

// Private constructor: instances used internally only
NavarroSeq::NavarroSeq(std::ifstream & ifs) 
{
	init(ifs);
}

void NavarroSeq::compress(string in_fname, string out_fname)
{

	std::ifstream ifs (in_fname);
	// Create private instance of NavarroSeq
	NavarroSeq* nseq = new NavarroSeq(ifs);

	string rfile = "temp_r_file";
	string ifile = "temp_i_file";
	string lfile = "temp_l_file";
	string nfile = "temp_n_file";

	string rmdr = nseq->parse_input_data(ifs, rfile, ifile, lfile, nfile);
	std::ofstream ofs(out_fname);
	nseq->print_header(ofs);

	std::ifstream ifs_r(rfile, std::ios_base::binary);
	std::ifstream ifs_i(ifile, std::ios_base::binary);
	std::ifstream ifs_l(lfile, std::ios_base::binary);
	std::ifstream ifs_n(nfile, std::ios_base::binary);

	ofs << ifs_r.rdbuf() << ifs_l.rdbuf() << ifs_i.rdbuf() << ifs_n.rdbuf();

	nseq->print_E_table_info(ofs);

	print_remainder(ofs);
	ofs << rmdr;

	remove("temp_r_file");
	remove("temp_i_file");
	remove("temp_l_file");
	remove("temp_n_file");
}

string NavarroSeq::decompress(string filename)
{
	
	std::ifstream ifs (filename, std::ifstream::in);

	unsigned int i,j,k;
	n = get_int(ifs, 32);
	r = get_int(ifs, 32);
	u = get_int(ifs, 32);
	num_blocks = get_int(ifs, 32);
	unsigned int num_combos = get_int(ifs, 32);

	// Skip stored alphabet
	for (i=0; i<r; i++) {
		char alpha = get_char(ifs, 8);
	}

	vector<unsigned int> r_vals (num_blocks, 0);
	vector<unsigned int> i_vals (num_blocks, 0);
	vector<unsigned int> l_partial_sums (num_blocks+1, 0);
	vector<unsigned int> n_partial_sums ((num_blocks+1)*r, 0);
	vector<unsigned int> E_table_depths (num_combos+1);
	vector<E_table*> E_table_ptrs; // Mapping from R index to E table

	// Populate r_vals
	unsigned int r_int_size = ceil(log2(num_combos));
	for (i=0; i<num_blocks; i++) {
		r_vals.at(i) = get_int(ifs, r_int_size);
	}

	flush_remainder();

	// Populate l_partial_sums
	large_block_size = ceil(log2(n*log2(r)));
	full_l_int_size = ceil(log2(n*log2(r)));
	relative_l_int_size = ceil(log2(u * ceil(log2(r)) * ceil(log2(n*log2(r)))));

	flush_remainder();
	unsigned int last_full_l_sum = 0;

	for (i=0; i<num_blocks+1; i++) {
		if (i % large_block_size == 0) {
			last_full_l_sum = get_int(ifs, full_l_int_size);
			l_partial_sums.at(i) = last_full_l_sum;
		} else {
			unsigned int relative_sum = get_int(ifs, relative_l_int_size);
			l_partial_sums.at(i) = relative_sum + last_full_l_sum;
		}
	}

	flush_remainder();

	// Populate i_vals
	for (i=0; i<num_blocks; i++) {
		if (l_partial_sums.at(i+1) == l_partial_sums.at(i)) {
			i_vals.at(i) = 0;
		}
		else {
			i_vals.at(i) = get_int(ifs, l_partial_sums.at(i+1)-l_partial_sums.at(i));
		}
	}

	flush_remainder();
	vector<unsigned int> last_full_n_sum(r,0);

	for (i=0; i<num_blocks+1; i++) {
		if (i % large_block_size == 0) {
			for(j=0; j<r; j++) {
				last_full_n_sum.at(j) = get_int(ifs, full_l_int_size);
				n_partial_sums.at(i*r + j) = last_full_n_sum.at(j);
			}
		} else {
			for(j=0; j<r; j++) {
				unsigned int relative_sum = get_int(ifs, relative_l_int_size);
				n_partial_sums.at(i*r+j) = relative_sum + last_full_n_sum.at(j);
			}
		}
	}

	flush_remainder();

	// Populate E Table Depths
	for (i=0; i<num_combos+1; i++) {
		E_table_depths.at(i) = get_int(ifs, 32);
	}

	flush_remainder();

	// Populate E Tables
	unsigned int e_table_rank_size = ceil(log2(u+1));
	for (i=0; i<num_combos; i++) {
		E_table* table = new E_table();
		unsigned int num_entries = E_table_depths.at(i+1) - E_table_depths.at(i);
		for (j=0; j<num_entries; j++) {
			char c_sequence[u+1];
			for (k=0; k<u; k++) {
				c_sequence[k] = get_char(ifs, 8);
			}
			c_sequence[u] = 0x00;
			string sequence (c_sequence);
			vector<char> ranks (r*u, 0);
			for (k=0; k<r*u; k++) {
				ranks.at(k) = get_char(ifs, e_table_rank_size);
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
	BitReader reader(&ifs);

	unsigned int u = reader.get_u();
	unsigned int block = floor(index/u);
	unsigned int l = index - block*u;
	unsigned int r_val = reader.get_r(block);
	unsigned int i_val = reader.get_i(block);
	unsigned char e_table_char = reader.get_E_Table_char(r_val, i_val, l);
	return e_table_char;
}

unsigned int NavarroSeq::rank(string fname, int index, char c)
{
	std::ifstream ifs (fname, std::ifstream::in);
	BitReader reader(&ifs);

	unsigned int u = reader.get_u();
	unsigned int block = floor(index/u);
	unsigned int l = index - block*u;
	unsigned int r_val = reader.get_r(block);
	unsigned int i_val = reader.get_i(block);
	unsigned int e_table_rank = reader.get_E_Table_rank(block, r_val, i_val, l, c);
	return e_table_rank;
}

unsigned int NavarroSeq::get_size(string fname) {
	std::ifstream ifs (fname, std::ifstream::in);
	return get_int(ifs, 32);
}
