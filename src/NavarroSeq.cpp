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
vector<unsigned int> curr_combination;
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
		//printf("Adding %02x to on_deck; on_deck now = %02x with %d bits.\n", c, on_deck, num_bits_on_deck);
	}
	else {
		int shift_length = num_bits_on_deck + num_bits - 8;
		//cout << "Shifting " << shift_length << " places." << endl;
		unsigned char shifted_c_filter = (unsigned char)0xff >> (8-(num_bits - shift_length));
		unsigned char shifted_c = (c >> shift_length) & shifted_c_filter;
		//printf("c = %02x; Shifted c = %02x using filter %02x\n", c, shifted_c, shifted_c_filter);
		unsigned char retval = on_deck + (shifted_c);
		//printf("retval = %02x + %02x = %02x\n", on_deck, shifted_c, retval);
		num_bits_on_deck = shift_length;
		on_deck = c << (8 - shift_length);
		ofs.put(retval);
		//printf("After printing, on_deck = %02x with %d bits.\n", on_deck, num_bits_on_deck);
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
		cout << "There aren't that many bits in an int." << endl;
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

	//printf("Int = %08x = %u\n", four_byte_union.integer, four_byte_union.integer);
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
		remaining_combo.at(i)--;
		for (j=i*u + prefix.length(); j<(i+1)*u; j++) {
			ranks.at(j)++;
		}
		get_etable_rows(next_prefix, ranks, remaining_combo, table);
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
			cout << "Alphabet " << percent_complete << "%% complete..." << endl;
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
	// Inefficient: we don't really care about compression efficiency, but this could be better
	// Figure out which index in R this is
	vector<unsigned int> current_combo (r,0); // Tranform string into character vector
	curr_combination = current_combo;
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
	//cout << "NUM BLOCKS = " << num_blocks << endl;
	print_int(ofs, combinations.size(), 32);
	/*for (int i=0; i<combinations.size(); i++) {
		cout << "Combinations(" << i << ") = <";
		vector<unsigned int> combo = combinations.at(i);
		for (int j=0; j<combo.size(); j++) {
			cout << combo.at(j);
			if (j != combo.size() - 1) cout << ",";
		}
		cout << ">" << endl;
	}
	*/
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

	//unsigned int total_e_table_size = 0;
	unsigned int e_table_rank_size = ceil(log2(u+1));
	for (size_t i=0; i<E_table_ptrs.size(); i++) {
		E_table* etable = E_table_ptrs.at(i);
		for (unsigned int j=0; j<etable->entries.size(); j++) {
			G_entry* entry = etable->entries.at(j);
			// Output compressed code
			for (unsigned int k=0; k<u; k++) {
				print_char(ofs, entry->sequence.at(k), 8);
				//total_e_table_size += 8;
			}
			for(unsigned int k=0; k<entry->ranks.size(); k++) {
				  //ofs.put(entry->ranks.at(k)); // log u bits
				print_char(ofs, entry->ranks.at(k), e_table_rank_size);
				//total_e_table_size += e_table_rank_size;
			}
		}
	}
}

string NavarroSeq::parse_input_data(std::ifstream & ifs, string r_fname, string i_fname, string l_fname, string n_fname) {
	
	unsigned int percent_complete=0, current_r=0, current_i=0, current_l=0, current_l_size=0, current_i_size=0, prev_l=0, last_full_l_sum = 0;
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

	for (size_t i=0; i<num_blocks; i++) { //Iterate over every full block
		if ((((float)i)/num_blocks)*100 > percent_complete + 1) {
			cout << "Compression " << percent_complete << "%% complete..." << endl;
			percent_complete++;
		}

		string current_block = get_next_block(ifs);
		current_r = get_r_index(current_block);
		E_table* table = E_table_ptrs.at(current_r);
		current_i = get_i_index(table, current_block);
		current_i_size = ceil(log2(table->entries.size()));
		
		//cout << "R(" << i << ") = " << current_r << " with " << r_int_size << " bits." << endl;
		r_printer.print_int(current_r, r_int_size);
		i_printer.print_int(current_i, current_i_size);

		/*
		* L and N partial sums are a bit trickier.  If this sequence begins a new L and N 'superblock', we store the full partial sums.
		* Otherwise, we store just the relative sum since the beginning of the current 'superblock'.
		*/
		if (i % large_block_size == 0) {
			current_l = (prev_l + current_i_size + last_full_l_sum);
			cout << "L(" << i << ") = " << current_l << " with " << full_l_int_size << " bits." << endl;
			l_printer.print_int(current_l, full_l_int_size);
			for (size_t j=0; j<r; j++) {
				unsigned int current_n = last_full_n_sum.at(j) + prev_n.at(j) + curr_combination.at(j);
				n_printer.print_int(current_n, 4);
				last_full_n_sum.at(j) = current_n;
				prev_n.at(j) = 0;
			}

			last_full_l_sum = current_l;
			prev_l = 0;
		} else {
			current_l = (prev_l + current_i_size);
			cout << "L(" << i << ") = " << current_l << " with " <<  relative_l_int_size << " bits." << endl;
			l_printer.print_int(current_l, relative_l_int_size);
			for (size_t j=0; j<r; j++) {
				unsigned int current_n = prev_n.at(j) + curr_combination.at(j);
				n_printer.print_int(current_n, 4);
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

	// Calculate bit sizes of each compressed element:
	r_int_size = ceil(log2(combinations.size()));
	large_block_size = ceil(log2(n*log2(r)));
	full_l_int_size = ceil(log2(n*log2(r)));
	relative_l_int_size = ceil(log2(u * ceil(log2(r)) * ceil(log2(n*log2(r)))));
	n_int_size = ceil(log2(n));
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

	cout << "Finished Successfully!" << endl;
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

	cout << "n = " << n << endl;
	cout << "r = " << r << endl;
	cout << "u = " << u << endl;
	cout << "num_blocks = " << num_blocks << endl;
	cout << "num_combos = " << num_combos << endl;

	// Skip stored alphabet
	cout << "Alphabet:" << endl;
	for (i=0; i<r; i++) {
		char alpha = get_char(ifs, 8);
		cout << "i = " << i << "; alpha = " << alpha << endl;
	}

	vector<unsigned int> r_vals (num_blocks, 0);
	vector<unsigned int> i_vals (num_blocks, 0);
	vector<unsigned int> l_partial_sums (num_blocks+1, 0);
	vector<unsigned int> E_table_depths (num_combos+1);
	vector<E_table*> E_table_ptrs; // Mapping from R index to E table

	// Populate r_vals
	cout << "R vals:" << endl;
	unsigned int r_int_size = ceil(log2(num_combos));
	for (i=0; i<num_blocks; i++) {
		r_vals.at(i) = get_int(ifs, r_int_size);
		cout << "i = " << i << "; R = " << r_vals.at(i) << endl;
	}

	/*
	for(i=0; i<10; i++) {
		unsigned char c = (unsigned char) ifs.get();
		printf("CHAR = %02x\n", c);
	}
	*/

	// Populate l_partial_sums
	cout << "L partial sums:" << endl;
	//unsigned int l_int_size = ceil(log2(n*log2(r)));
	large_block_size = ceil(log2(n*log2(r)));
	full_l_int_size = ceil(log2(n*log2(r)));
	relative_l_int_size = ceil(log2(n*log2(r)));

	unsigned int last_full_l_sum = 0;

	for (i=0; i<num_blocks+1; i++) {
		if (i % large_block_size == 0) {
			last_full_l_sum = get_int(ifs, full_l_int_size);
			l_partial_sums.at(i) = last_full_l_sum;
			cout << "L(" << i << ") = " << l_partial_sums.at(i) << "[" << full_l_int_size << "] bits." << endl;
		} else {
			unsigned int relative_sum = get_int(ifs, relative_l_int_size);
			l_partial_sums.at(i) = relative_sum + last_full_l_sum;
			cout << "L(" << i << ") = " << relative_sum << " + " << last_full_l_sum << " = " << l_partial_sums.at(i) << "[" << relative_l_int_size << "] bits." << endl;
		}
	}

	/*
	// Populate i_vals
	cout << "I vals:" << endl;
	for (i=0; i<num_blocks; i++) {
		if (l_partial_sums.at(i+1) == l_partial_sums.at(i)) i_vals.at(i) = 0;
		else i_vals.at(i) = get_int(ifs, l_partial_sums.at(i+1)-l_partial_sums.at(i));
		cout << "  " << i_vals.at(i) << endl;
	}

	// Skip n_partial_sums
	//ifs.ignore((num_blocks+1)*r*4+4, EOF); 
	unsigned int n_int_size = ceil(log2(n));
	for (i=0; i<(num_blocks+1)*r; i++) {
		get_int(ifs, n_int_size);
	}

	// Populate E Table Depths
	for (i=0; i<num_combos+1; i++) {
		E_table_depths.at(i) = get_int(ifs, 32);
	}

	// Populate E Tables
	unsigned int e_table_rank_size = ceil(log2(u+1));
	for (i=0; i<num_combos; i++) {
		E_table* table = new E_table();
		unsigned int num_entries = E_table_depths.at(i+1) - E_table_depths.at(i);
		for (j=0; j<num_entries; j++) {
			char c_sequence[u+1];
			cout << "String: " << endl;
			for (k=0; k<u; k++) {
				c_sequence[k] = get_char(ifs, 8);
				cout << "  char = " << c_sequence[k] << endl;
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
		cout << "S = " << entry->sequence << endl;
	}

	cout << "DECOMPRESSED STRING = " << s << rmdr << endl;
	
	return s;
	*/
	return "Hai!";
}

char NavarroSeq::access(string fname, int index)
{
	/*
	std::ifstream ifs (fname, std::ifstream::in);

	unsigned int i,j,k;
	unsigned int n = get_int(ifs, 32);
	unsigned int r = get_int(ifs, 32);
	unsigned int u = get_int(ifs, 32);
	unsigned int num_blocks = get_int(ifs, 32);
	unsigned int num_combos = get_int(ifs, 32);

	// Skip stored alphabet
	for (i=0; i<r; i++) {
		//cout << i << endl;
		ifs.get();
	}

	unsigned int block = floor(index/u);
	unsigned int l = index - block*u;

	ifs.seekg((5 + block) * 4 + r);
	unsigned int r_val = get_int(ifs, 32);

	ifs.seekg((5 + num_blocks + block) * 4 + r);
	unsigned int i_val = get_int(ifs, 32);

	ifs.seekg((7+num_blocks*(3+r) + r + r_val)*4 + r);
	unsigned int E_depth = get_int(ifs, 32);

	ifs.seekg((8+num_blocks*(3+r)+r+num_combos)*4 + E_depth*((r + 1) * u) + l + r);
	char retval = ifs.get();

	return retval;
	*/
	return 0x00;
}

unsigned int NavarroSeq::rank(string fname, char c, int index)
{
	/*
	std::ifstream ifs (fname, std::ifstream::in);

	unsigned int char_index_in_alphabet = 0; // hard coded for now

	unsigned int i,j,k;
	unsigned int n = get_int(ifs, 32);
	unsigned int r = get_int(ifs, 32);
	unsigned int u = get_int(ifs, 32);
	unsigned int num_blocks = get_int(ifs, 32);
	unsigned int num_combos = get_int(ifs, 32);

	// Find in stored alphabet
	for (i=0; i<r; i++) {
		char current_char = ifs.get();
		if (current_char == c) char_index_in_alphabet = i;
	}

	unsigned int block = floor(index/u);
	unsigned int l = index - block*u;

	ifs.seekg((5 + block) * 4 + r);
	unsigned int r_val = get_int(ifs, 32);

	ifs.seekg((5 + num_blocks + block) * 4 + r);
	unsigned int i_val = get_int(ifs, 32);

	ifs.seekg((7 + 3*num_blocks + r*block + char_index_in_alphabet) * 4 + r);
	unsigned int n_partial_sum = get_int(ifs, 32);

	ifs.seekg((7+num_blocks*(3+r)+r + r_val)*4 + r);
	unsigned int E_depth = get_int(ifs, 32);

	ifs.seekg((8+num_blocks*(3+r)+r+num_combos)*4 + E_depth*((r + 1) * u) + u + u*char_index_in_alphabet + l + r);
	char block_rank = ifs.get();

	unsigned int rank = n_partial_sum + block_rank;

	//cout << "RANK for index " << index << " = " << n_partial_sum << " + " << (int)block_rank << " = " << rank << endl;

	return rank;
	*/
	return 4;
}

unsigned int NavarroSeq::select(string fname, char c, int index)
{
	return 55;
}

unsigned int NavarroSeq::get_size(string fname) {
	std::ifstream ifs (fname, std::ifstream::in);
	return get_int(ifs, 32);
}

int main(int argc, char** argv) {

	bool compress = false;
	bool decompress = false;

	string fileIn = "";
	string fileOut = "";

	if (argc > 1) {
		if (strcmp("-compress", argv[1]) == 0) {
			compress = true;
			fileIn = argv[2];
			fileOut = argv[3];
		} 
		else if (strcmp("-decompress", argv[1]) == 0) {
			decompress = true;
			fileIn = argv[2];
		} else {
			return 0;
		}
	}

	if (compress) {
		NavarroSeq::compress(fileIn,  fileOut);
	}
	else if (decompress) {
		NavarroSeq::decompress(fileIn);
	}
/*
	unsigned int testint = 42;
	std::ofstream ofs("test-byte-compression.txt");
	print_char_with_num_bits(ofs, 0x03, 2);
	print_char_with_num_bits(ofs, 0x00, 2);
	print_char_with_num_bits(ofs, 0x02, 2);
	print_char_with_num_bits(ofs, 0x6f, 6);
	print_with_num_bits(ofs, 0xaabbccdd, 32);
	print_with_num_bits(ofs, testint, 7);
	print_with_num_bits(ofs, 0x1abbccdd, 29);
	print_remainder(ofs);
*/

/*
	std::ifstream ifs ("test-byte-compression.txt", std::ifstream::in);
	get_char_with_num_bits(ifs, 2);
	get_char_with_num_bits(ifs, 2);
	get_char_with_num_bits(ifs, 2);
	get_char_with_num_bits(ifs, 6);
	get_int_with_num_bits(ifs, 32);
	get_int_with_num_bits(ifs, 7);
	get_int_with_num_bits(ifs, 29);
*/

	return 0;
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
