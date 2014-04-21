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
	unsigned long get_long(int num_bits);
	void flush_remainder();
	unsigned int get_u();
	unsigned int get_r(unsigned long block);
	unsigned long get_l(unsigned long block);
	unsigned int get_i(unsigned long block);
	unsigned long get_n(unsigned long block, char c);
	unsigned char get_E_Table_char(unsigned int r_val, unsigned int i_val, unsigned int index);
	unsigned long get_E_Table_rank(unsigned long block, unsigned int r_val, unsigned int i_val, unsigned int index, char c);
private:
	unsigned int get_stored_int(unsigned long header_bytes, unsigned long preceding_bits, unsigned int read_size);
	unsigned long get_stored_long(unsigned long header_bytes, unsigned long preceding_bits, unsigned int read_size);
	unsigned int get_char_index_in_alphabet(char c);
	unsigned long get_r_bytes();
	unsigned long get_l_bytes();
	unsigned long get_i_bytes();
	unsigned long get_n_bytes();
	unsigned long get_E_Table_depth(unsigned int r_val);
	std::ifstream* ifs;
	unsigned char on_deck;
	unsigned int num_bits_on_deck;
	unsigned long n;
	unsigned int r;
	unsigned int u;
	unsigned long num_blocks;
	unsigned int num_combos;
	unsigned int r_size;
	unsigned int superblock_size;
	unsigned int full_partial_sum_size;
	unsigned int relative_partial_sum_size;
	unsigned int e_table_rank_size;
	unsigned long total_bytes_through_r;
	unsigned long total_bytes_through_l;
	unsigned long total_bytes_through_i;
	unsigned long total_bytes_through_n;
	union {
		unsigned int integer;
		char byte[4];
	} four_byte_union;
	union {
		unsigned int long_val;
		char byte[8];
	} eight_byte_union;
};

#endif