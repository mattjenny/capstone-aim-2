#include "BitReader.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>

using namespace std;

BitReader::BitReader(std::ifstream * ifs_in) {
	ifs = ifs_in;
	on_deck = 0x00;
	num_bits_on_deck = 0;
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