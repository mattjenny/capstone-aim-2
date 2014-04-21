#ifndef BITPRINTER_H
#define BITPRINTER_H

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

class BitPrinter
{
public:
	BitPrinter(std::ofstream * ofs_in);
	void print_char(char c_in, unsigned int num_bits);
	void print_int(unsigned long integer, unsigned int num_bits);
	void print_remainder();
private:
	std::ofstream* ofs;
	char on_deck;
	unsigned int num_bits_on_deck;
	union {
		unsigned long long_val;
		char byte[8];
	} eight_byte_union;
};

#endif