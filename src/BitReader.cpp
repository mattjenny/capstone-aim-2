#include "BitReader.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <math.h>

using namespace std;

unsigned int BYTES_IN_INT = 4;
unsigned int HEADER_SIZE = 5*BYTES_IN_INT;

BitReader::BitReader(std::ifstream * ifs_in) {
	ifs = ifs_in;
	on_deck = 0x00;
	num_bits_on_deck = 0;
	n = get_int(32);
	r = get_int(32);
	u = get_int(32);
	num_blocks = get_int(32);
	num_combos = get_int(32);

	r_size = ceil(log2(num_combos));
	superblock_size = ceil(log2(n*log2(r)));
	full_partial_sum_size = ceil(log2(n*log2(r)));
	relative_partial_sum_size = ceil(log2(u * ceil(log2(r)) * ceil(log2(n*log2(r)))));
	e_table_rank_size = ceil(log2(u+1));

	total_bytes_through_r = get_r_bytes();
	total_bytes_through_l = get_l_bytes();
	total_bytes_through_i = get_i_bytes();
	total_bytes_through_n = get_n_bytes();

}

unsigned char BitReader::get_char(int num_bits) {
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
		unsigned char c = (unsigned char) ifs->get();
		retval = (on_deck >> (8 - num_bits)) + (c >> (8 - num_bits + num_bits_on_deck));
		on_deck = c << (num_bits - num_bits_on_deck);
		num_bits_on_deck = 8 - num_bits + num_bits_on_deck;
	}

	return retval;
}

unsigned int BitReader::get_int(int num_bits) {
	if (num_bits > 32) {
		cout << "There aren't " << num_bits << " bits in an int." << endl;
		return 0;
	} else if (num_bits == 0) {
		return 0;
	}

	for (int i=0; i<4; i++) {
		if ((3-i)*8 >= num_bits) {
			four_byte_union.byte[3-i] = 0x00;
		}
		else if ((4-i)*8 <= num_bits) {
			four_byte_union.byte[3-i] = get_char(8);
		} 
		else {
			four_byte_union.byte[3-i] = get_char(num_bits % 8);
		}
	}

	return four_byte_union.integer;
}

void BitReader::flush_remainder() {
	on_deck = 0x00;
	num_bits_on_deck = 0;
}

unsigned int BitReader::get_stored_int(unsigned int header_bytes, unsigned int preceding_bits, unsigned int read_size) 
{
	flush_remainder();

	unsigned int preceding_bytes = preceding_bits/8;
	unsigned int preceding_rmdr = preceding_bits % 8;

	ifs->seekg(header_bytes + preceding_bytes);
	get_int(preceding_rmdr);
	unsigned int retval = get_int(read_size);
	return retval;
}

unsigned int BitReader::get_r_bytes() {
	return HEADER_SIZE + r + ceil(num_blocks*r_size/8.0);
}

unsigned int BitReader::get_l_bytes() {
	unsigned int num_r_bytes = total_bytes_through_r;
	unsigned int num_superblocks = ((num_blocks) / superblock_size) + 1;
	unsigned int num_relative_blocks = num_blocks - num_superblocks + 1;
	unsigned int num_l_bytes = ceil((num_superblocks*full_partial_sum_size + num_relative_blocks*relative_partial_sum_size)/8.0);
	return num_r_bytes + num_l_bytes;
}

unsigned int BitReader::get_i_bytes() {
	unsigned int num_l_bytes = total_bytes_through_l;
	unsigned int num_i_bits = get_l(num_blocks);
	unsigned int num_i_bytes = ceil(num_i_bits/8.0);
	return num_l_bytes + num_i_bytes;
}

unsigned int BitReader::get_n_bytes() {
	unsigned int num_i_bytes = total_bytes_through_i;
	unsigned int num_r_bytes = total_bytes_through_r;
	unsigned int num_superblocks = ((num_blocks) / superblock_size) + 1;
	unsigned int num_relative_blocks = num_blocks - num_superblocks + 1;
	unsigned int num_n_bytes = ceil((num_superblocks*full_partial_sum_size + num_relative_blocks*relative_partial_sum_size)*r/8.0);
	return num_i_bytes + num_n_bytes;
}

unsigned int BitReader::get_char_index_in_alphabet(char c) {
	ifs->seekg(0);
	for (size_t i=0; i<HEADER_SIZE; i++) ifs->get(); //skip first 5 ints
	for (size_t i=0; i<r; i++) {
		char current_char = ifs->get();
		if (current_char == c) return i;
	}
	return 42; //should add error checking
}

unsigned int BitReader::get_u()
{
	return u;
}

unsigned int BitReader::get_r(unsigned int block) 
{
	return get_stored_int(HEADER_SIZE + r, block*r_size, r_size);
}

unsigned int BitReader::get_l(unsigned int block) 
{
	unsigned int num_r_bytes = total_bytes_through_r;
	unsigned int num_bits_in_superblock = full_partial_sum_size + (superblock_size-1)*relative_partial_sum_size;
	unsigned int superblock = block / superblock_size;

	unsigned int full_partial_sum = get_stored_int(num_r_bytes, superblock*num_bits_in_superblock, full_partial_sum_size);
	unsigned int relative_partial_sum;
	if (block % superblock_size == 0)  relative_partial_sum = 0;
	else {
		unsigned int preceding_superblocks = ((block - 1) / superblock_size) + 1;
		unsigned int preceding_relative_blocks = block - preceding_superblocks;
		unsigned int preceding_bits = preceding_superblocks*full_partial_sum_size + preceding_relative_blocks*relative_partial_sum_size;
		relative_partial_sum = get_stored_int(num_r_bytes, preceding_bits, relative_partial_sum_size);
	}

	return full_partial_sum + relative_partial_sum;
}

unsigned int BitReader::get_i(unsigned int block) 
{
	unsigned int num_l_bytes = total_bytes_through_l;
	unsigned int prev_i_bits = get_l(block);
	unsigned int current_i_size = get_l(block+1) - prev_i_bits;
	unsigned int retval = get_stored_int(num_l_bytes, prev_i_bits, current_i_size);
	return retval;
}

unsigned int BitReader::get_n(unsigned int block, char c) 
{
	unsigned int previous_bytes = total_bytes_through_i;
	unsigned int char_index = get_char_index_in_alphabet(c);

	unsigned int num_bits_in_superblock = (full_partial_sum_size + (superblock_size-1)*relative_partial_sum_size)*r;
	unsigned int superblock = block / superblock_size;
	unsigned int previous_bits = (superblock*num_bits_in_superblock)+(full_partial_sum_size*char_index);
	unsigned int full_partial_sum = get_stored_int(previous_bytes, previous_bits, full_partial_sum_size);

	unsigned int relative_partial_sum;
	if (block % superblock_size == 0)  relative_partial_sum = 0;
	else {
		unsigned int preceding_superblocks = ((block - 1) / superblock_size) + 1;
		unsigned int preceding_relative_blocks = block - preceding_superblocks;
		unsigned int preceding_bits = (preceding_superblocks*full_partial_sum_size + preceding_relative_blocks*relative_partial_sum_size)*r + char_index*relative_partial_sum_size;
		relative_partial_sum = get_stored_int(previous_bytes, preceding_bits, relative_partial_sum_size);
	}
	return full_partial_sum + relative_partial_sum;
}

unsigned int BitReader::get_E_Table_depth(unsigned int r_val)
{
	unsigned int previous_bytes = total_bytes_through_n;
	return get_stored_int(previous_bytes, r_val*32, 32);
}

unsigned char BitReader::get_E_Table_char(unsigned int r_val, unsigned int i_val, unsigned int index) 
{
	unsigned int E_Table_depth = get_E_Table_depth(r_val);
	unsigned int previous_bytes = total_bytes_through_n + 4*(num_combos + 1);// WHY 3??  I don't know.  It's magic.
	unsigned int previous_bits = (E_Table_depth + i_val)*((r * u)*e_table_rank_size + u*8) + 8*index;
	unsigned char retval = (unsigned char) get_stored_int(previous_bytes, previous_bits, 8);
	return retval;
}

unsigned int BitReader::get_E_Table_rank(unsigned int block, unsigned int r_val, unsigned int i_val, unsigned int index, char c) 
{
	unsigned int char_index = get_char_index_in_alphabet(c);
	unsigned int E_Table_depth = get_E_Table_depth(r_val);
	unsigned int previous_bytes = total_bytes_through_n + 4*(num_combos + 1);
	unsigned int previous_bits = (E_Table_depth + i_val)*((r * u)*e_table_rank_size + u*8) + 8*u + (char_index*u + index)*e_table_rank_size;
	unsigned int block_rank = get_stored_int(previous_bytes, previous_bits, e_table_rank_size);
	unsigned int n = get_n(block, c);
	return n + block_rank;
}