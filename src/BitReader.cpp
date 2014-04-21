#include "BitReader.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <math.h>

using namespace std;

unsigned int BYTES_IN_INT = 4;
unsigned int HEADER_SIZE = 7*BYTES_IN_INT;

BitReader::BitReader(std::ifstream * ifs_in) {
	ifs = ifs_in;
	on_deck = 0x00;
	num_bits_on_deck = 0;
	//ifs->seekg(0, ifs->end);
	n = get_long(64);
	r = get_int(32);
	u = get_int(32);
	num_blocks = get_long(64);
	num_combos = get_int(32);

	cout << "n = " << n << endl;
	cout << "r = " << r << endl;
	cout << "alphabet = " << endl;
	ifs->seekg(0);
	for (size_t i=0; i<HEADER_SIZE; i++) ifs->get(); //skip first 5 ints
	for (size_t i=0; i<r; i++) {
		char current_char = ifs->get();
		printf("   %02x\n", current_char);
	}
	cout << "u = " << u << endl;
	cout << "num blocks = " << num_blocks << endl;
	cout << "num_combos = " << num_combos << endl;

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

unsigned long BitReader::get_long(int num_bits) {
	if (num_bits > 64) {
		cout << "There aren't " << num_bits << " bits in a long." << endl;
		return 0;
	} else if (num_bits == 0) {
		return 0;
	}

	for (int i=0; i<8; i++) {
		if ((7-i)*8 >= num_bits) {
			eight_byte_union.byte[7-i] = 0x00;
		}
		else if ((8-i)*8 <= num_bits) {
			eight_byte_union.byte[7-i] = get_char(8);
		} 
		else {
			eight_byte_union.byte[7-i] = get_char(num_bits % 8);
		}
	}

	return eight_byte_union.long_val;
}

void BitReader::flush_remainder() {
	on_deck = 0x00;
	num_bits_on_deck = 0;
}

unsigned int BitReader::get_stored_int(unsigned long header_bytes, unsigned long preceding_bits, unsigned int read_size) 
{
	cout << "Calling get_stored_int with " << header_bytes << " header bytes, " << preceding_bits << " preceding bits, and " << read_size << " bits read size." << endl;
	flush_remainder();

	unsigned long preceding_bytes = preceding_bits/8;
	unsigned long preceding_rmdr = preceding_bits % 8;

	ifs->seekg(header_bytes + preceding_bytes);
	get_int(preceding_rmdr);
	unsigned int retval = get_int(read_size);
	cout << "preceding bytes = " << preceding_bytes << "; rmdr = " << preceding_rmdr << endl;
	cout << "  retval = " << retval << endl;
	return retval;
}

unsigned long BitReader::get_stored_long(unsigned long header_bytes, unsigned long preceding_bits, unsigned int read_size) 
{
	cout << "Calling get_stored_int with " << header_bytes << " header bytes, " << preceding_bits << " preceding bits, and " << read_size << " bits read size." << endl;
	flush_remainder();

	unsigned long preceding_bytes = preceding_bits/8;
	unsigned long preceding_rmdr = preceding_bits % 8;

	ifs->seekg(header_bytes + preceding_bytes);
	get_int(preceding_rmdr);
	unsigned long retval = get_long(read_size);
	cout << "preceding bytes = " << preceding_bytes << "; rmdr = " << preceding_rmdr << endl;
	cout << "  retval = " << retval << endl;
	return retval;
}

unsigned long BitReader::get_r_bytes() {
	cout << "GETTING R BYTES"  << endl;
	return HEADER_SIZE + r + ceil(num_blocks*r_size/8.0);
}

unsigned long BitReader::get_l_bytes() {
	cout << "GETTING L BYTES" << endl;
	unsigned int num_r_bytes = total_bytes_through_r;
	unsigned int num_superblocks = ((num_blocks - 1) / superblock_size) + 1;
	cout << "Num superblocks = " << num_superblocks << endl;
	unsigned int num_relative_blocks = num_blocks - num_superblocks;
	cout << "Num relative blocks = " << num_relative_blocks << endl;
	unsigned int num_l_bytes = ceil((num_superblocks*full_partial_sum_size + num_relative_blocks*relative_partial_sum_size)/8.0);
	cout << "L bytes = " << num_r_bytes << " + " << num_l_bytes << " = " << num_r_bytes + num_l_bytes << endl;
	return num_r_bytes + num_l_bytes;
}

unsigned long BitReader::get_i_bytes() {
	cout << "GETTING I BYTES" << endl;
	unsigned int num_l_bytes = total_bytes_through_l;
	unsigned int num_i_bits = get_l(num_blocks);
	unsigned int num_i_bytes = ceil(num_i_bits/8.0);
	return num_l_bytes + num_i_bytes;
}

unsigned long BitReader::get_n_bytes() {
	cout << "GETTING N BYTES" << endl;
	unsigned int num_i_bytes = total_bytes_through_i;
	unsigned int num_r_bytes = total_bytes_through_r;
	unsigned int num_superblocks = ((num_blocks - 1) / superblock_size) + 1;
	unsigned int num_relative_blocks = num_blocks - num_superblocks;
	unsigned int num_n_bytes = ceil((num_superblocks*full_partial_sum_size + num_relative_blocks*relative_partial_sum_size)*r/8.0);
	cout << "num i bytes = " << num_i_bytes << "; num n bytes = " << num_n_bytes << endl;
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

unsigned int BitReader::get_r(unsigned long block) 
{
	cout << "GETTING R" << endl;
	return get_stored_int(HEADER_SIZE + r, block*r_size, r_size);
}

unsigned long BitReader::get_l(unsigned long block) 
{
	cout << "GETTING L" << endl;
	unsigned long num_r_bytes = total_bytes_through_r;
	unsigned int num_bits_in_superblock = full_partial_sum_size + (superblock_size-1)*relative_partial_sum_size;
	unsigned long superblock = block / superblock_size;

	cout << "Calculating L full partial sum:" << endl;
	unsigned long full_partial_sum = get_stored_long(num_r_bytes, superblock*num_bits_in_superblock, full_partial_sum_size);
	unsigned int relative_partial_sum;
	if (block % superblock_size == 0)  relative_partial_sum = 0;
	else {
		unsigned int preceding_superblocks = ((block - 1) / superblock_size) + 1;
		unsigned int preceding_relative_blocks = block - preceding_superblocks;
		unsigned int preceding_bits = preceding_superblocks*full_partial_sum_size + preceding_relative_blocks*relative_partial_sum_size;
		cout << "Calculating L relative partial sum:" << endl;
		relative_partial_sum = get_stored_int(num_r_bytes, preceding_bits, relative_partial_sum_size);
	}

	return full_partial_sum + relative_partial_sum;
}

unsigned int BitReader::get_i(unsigned long block) 
{
	cout << "GETTING I" << endl;
	unsigned long num_l_bytes = total_bytes_through_l + 1;

	unsigned long prev_i_bits = get_l(block);
	unsigned long current_i_size = get_l(block+1) - prev_i_bits;
	ifs->seekg(num_l_bytes);
	for(int i=0; i<20; i++) {
		printf("+++++ %02x\n", ifs->get());
	}
	cout << "Getting I..." << endl;
	unsigned int retval = get_stored_int(num_l_bytes, prev_i_bits, current_i_size);
	cout << "Retval = " << retval << endl;
	return retval;
}

unsigned long BitReader::get_n(unsigned long block, char c) 
{
	cout << "GETTING N" << endl;
	unsigned long previous_bytes = total_bytes_through_i + 1;
	unsigned int char_index = get_char_index_in_alphabet(c);

	unsigned int num_bits_in_superblock = (full_partial_sum_size + (superblock_size-1)*relative_partial_sum_size)*r;
	unsigned long superblock = block / superblock_size;
	unsigned long previous_bits = (superblock*num_bits_in_superblock)+(full_partial_sum_size*char_index);
	unsigned long full_partial_sum = get_stored_long(previous_bytes, previous_bits, full_partial_sum_size);

	unsigned int relative_partial_sum;
	if (block % superblock_size == 0)  relative_partial_sum = 0;
	else {
		unsigned int preceding_superblocks = ((block - 1) / superblock_size) + 1;
		unsigned int preceding_relative_blocks = block - preceding_superblocks;
		unsigned int preceding_bits = (preceding_superblocks*full_partial_sum_size + preceding_relative_blocks*relative_partial_sum_size)*r + char_index*relative_partial_sum_size;
		relative_partial_sum = get_stored_int(previous_bytes, preceding_bits, relative_partial_sum_size);
	}

	cout << " ----> N full = " << full_partial_sum << " and partial = " << relative_partial_sum << endl;
	return full_partial_sum + relative_partial_sum;
}

unsigned long BitReader::get_E_Table_depth(unsigned int r_val)
{
	cout << "GETTING E TABLE DEPTH with r_val = " << r_val << endl;
	unsigned long previous_bytes = total_bytes_through_n;
	cout << "  previous bytes = " << previous_bytes << endl;
	return get_stored_long(previous_bytes+8, r_val*32, 64);
}

unsigned char BitReader::get_E_Table_char(unsigned int r_val, unsigned int i_val, unsigned int index) 
{
	unsigned int E_Table_depth = get_E_Table_depth(r_val);
	unsigned int previous_bytes = total_bytes_through_n + 4*(num_combos+2);
	unsigned int previous_bits = (E_Table_depth + i_val)*((r * u)*e_table_rank_size + u*8) + 8*index;
	ifs->seekg(total_bytes_through_n + 8 + 4*num_combos);
	unsigned int prevbits = 0;
	for (int i=0; i<10; i++) {
		for (int j=0; j<4; j++) {
			char c = get_char(8);
			if(prevbits == previous_bits) printf("-------> ");
			printf("!!! C = %02x\n", c);
			prevbits += 8;
		}
		for (int j=0; j<r*u; j++) {
			char c = get_char(e_table_rank_size);
			if(prevbits == previous_bits) printf("-------> ");
			printf("^^^ C = %02x\n", c);
			prevbits += e_table_rank_size;
		}
	}
	unsigned char retval = (unsigned char) get_stored_int(previous_bytes, previous_bits, 8);
	printf("Retrieved E Table Char: %02x\n", retval);
	return retval;
}

unsigned long BitReader::get_E_Table_rank(unsigned long block, unsigned int r_val, unsigned int i_val, unsigned int index, char c) 
{
	unsigned int char_index = get_char_index_in_alphabet(c);
	unsigned long E_Table_depth = get_E_Table_depth(r_val);
	unsigned long previous_bytes = total_bytes_through_n + 8*(num_combos+2);
	unsigned long previous_bits = (E_Table_depth + i_val)*((r * u)*e_table_rank_size + u*8) + 8*u + (char_index*u + index)*e_table_rank_size;
	unsigned int block_rank = get_stored_int(previous_bytes, previous_bits, e_table_rank_size);
	printf("Retrieved E Table Rank: %u\n", block_rank);
	unsigned long n = get_n(block, c);
	cout << "N = " << n << endl;
	return n + block_rank;
}