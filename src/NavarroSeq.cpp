#include "NavarroSeq.h"
#include <string>
#include <math.h>
#include <unordered_set>
#include <stdio.h>
#include <vector>

using namespace std;

string NavarroSeq::compress(string s)
{
	size_t n = s.size();
	size_t i,r;

	unordered_set<char> unique_chars;
	for (i=0; i<n; i++) {
		unique_chars.insert(s.at(i));
	}
	r = unique_chars.size();

	size_t u = floor(0.5*log((double) n)/log((double) r));
	printf("n = %i; r = %i; u = %i\n", (int)n, (int)r, (int)u);

	//example enumeration (while it's still fresh on my mind) :

	vector<vector<unsigned int> > combinations = get_all_combinations(unique_chars.size(), u);
	printf ("Vector size is %u\n", (unsigned int) combinations.size());

	// get a chunk of u characters
	// make a vector of length r containing numbers of each character

	// I
	// R
	// E tables
	// F mapping

	return "Hello World!";
}

string NavarroSeq::decompress(string s)
{
	return "Hello again, world!";
}

char NavarroSeq::access(string s, int index)
{
	return 4;
}

unsigned int NavarroSeq::rank(string s, char c, int index)
{
	return 14;
}

unsigned int NavarroSeq::select(string s, char c, int index)
{
	return 55;
}

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
			values.at(i)++;
			enumerate(values, combos, i, r, u-1);
		}
	}
}

int main() {
	// test sequence
	NavarroSeq::compress("abbabaababbababcabbabaababbababcabbabaababbababcabbabaababbababcabbabaababbababcabc");
	return 0;
}