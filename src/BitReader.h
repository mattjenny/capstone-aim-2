#ifndef BITREADER_H
#define BITREADER_H

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

class BitReader
{
public:
	BitReader(std::ifstream * ifs_in);
	unsigned char get_char(int num_bits);
	unsigned int get_int(int num_bits);	
	void flush_remainder();
	unsigned int get_u();
	unsigned int get_r(unsigned int block);
	unsigned int get_l(unsigned int block);
	unsigned int get_i(unsigned int block);
	unsigned int get_n(unsigned int block, char c);
	unsigned char get_E_Table_char(unsigned int r_val, unsigned int i_val, unsigned int index);
	unsigned int get_E_Table_rank(unsigned int block, unsigned int r_val, unsigned int i_val, unsigned int index, char c);
private:
	unsigned int get_stored_int(unsigned int header_bytes, unsigned int preceding_bits, unsigned int read_size);
	unsigned int get_char_index_in_alphabet(char c);
	unsigned int get_r_bytes();
	unsigned int get_l_bytes();
	unsigned int get_i_bytes();
	unsigned int get_n_bytes();
	unsigned int get_E_Table_depth(unsigned int r_val);
	std::ifstream* ifs;
	unsigned char on_deck;
	unsigned int num_bits_on_deck;
	unsigned int n;
	unsigned int r;
	unsigned int u;
	unsigned int num_blocks;
	unsigned int num_combos;
	unsigned int r_size;
	unsigned int superblock_size;
	unsigned int full_partial_sum_size;
	unsigned int relative_partial_sum_size;
	unsigned int e_table_rank_size;
	unsigned int total_bytes_through_r;
	unsigned int total_bytes_through_l;
	unsigned int total_bytes_through_i;
	unsigned int total_bytes_through_n;
	union {
		unsigned int integer;
		char byte[4];
	} four_byte_union;
};

#endif