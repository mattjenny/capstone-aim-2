#include "BitPrinter.h"
#include <stdio.h>

using namespace std;

BitPrinter::BitPrinter(std::ofstream* ofs_in) {
	ofs = ofs_in;
	on_deck = 0x00;
	num_bits_on_deck = 0;
}

void BitPrinter::print_char(char c_in, unsigned int num_bits) {
	unsigned char c = (unsigned char) c_in;
	if (num_bits == 0) return;
	if (num_bits_on_deck + num_bits < 8) {
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
		ofs->put(retval);
		//printf("After printing, on_deck = %02x with %d bits.\n", on_deck, num_bits_on_deck);
	}
}

void BitPrinter::print_int(unsigned int integer, unsigned int num_bits) {
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
	print_char(four_byte_union.byte[max_byte], mod_bits);

	for (i=max_byte - 1; i>=0; i--) {
		print_char(four_byte_union.byte[i], 8);
	} 
}

int main() {
	return 0;
}