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
private:
	std::ifstream* ifs;
	char on_deck;
	unsigned int num_bits_on_deck;
	union {
		unsigned int integer;
		char byte[4];
	} four_byte_union;
};

#endif