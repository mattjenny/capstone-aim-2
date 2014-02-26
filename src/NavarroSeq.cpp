#include "NavarroSeq.h"
#include <string>
#include <math.h>
#include <unordered_set>
#include <stdio.h>
#include <vector>
#include <list>
#include <limits>
#include <iostream>
#include <sstream>

using namespace std;

size_t nseq_n, nseq_r, nseq_u;
vector<char> nseq_alphabet;

// Private constructor: instances used internally only
NavarroSeq::NavarroSeq(size_t n, size_t r, size_t u, list<char> alphabet) 
{
	nseq_n = n;
	nseq_r = r;
	nseq_u = u;
	for (list<char>::iterator it = alphabet.begin(); it != alphabet.end(); ++it) {
		nseq_alphabet.push_back(*it);
	}
}

string NavarroSeq::toBinStr(unsigned int val)
{
  string binStr;
  unsigned int mask = 1 << (sizeof(int) * 8 - 1);

   for(int i = 0; i < sizeof(int) * 8; i++)
   {
      if( (val & mask) == 0 )
         binStr += '0' ;
      else
         binStr += '1' ;

      mask  >>= 1;
   }

   std::size_t pos =  binStr.find_first_not_of("0");
   if(pos > 0)
     binStr.erase(0,pos); 

   //cout << val << " " << binStr << endl;
   return binStr;
}

string NavarroSeq::compress(string s)
{
	/*
	* n is the length of the full sequence s
	* r is the number of characters in the alphabet
	* u is the block size
	* num_blocks is the number of u-length blocks in the sequence
	* i,j,k are counters used throughout this method
	*/
	size_t n,r,u,num_blocks,i,j,k;
	n = s.size();

	// Figure out the alphabet that this sequence uses
	unordered_set<char> unique_chars;
	for (i=0; i<n; i++) {
		unique_chars.insert(s.at(i));
	}
	r = unique_chars.size();

	// Sort the alphabet
	list<char> alphabet;
	for(unordered_set<char>::iterator set_it = unique_chars.begin(); set_it != unique_chars.end(); ++set_it) {
		alphabet.push_back(*set_it);
	}
	alphabet.sort();

	// Calculate block length u and number of blocks
	u = floor(0.5*log((double) n)/log((double) r));
	num_blocks = ceil((n+0.0)/u);

	// Create private instance of NavarroSeq
	NavarroSeq* nseq = new NavarroSeq(n, r, u, alphabet);
	
	// Enumerate and index all possible COMBINATIONS of alphabet characters
	vector<vector<unsigned int> > combinations = get_all_combinations(unique_chars.size(), u);

	//Create E tables:
	vector<E_table*> E_table_ptrs; // Mapping from R index to E table
	for (i=0; i<combinations.size(); i++) {
		// Create E table for R[i]
		E_table* etable = nseq->get_etable(combinations.at(i));
		E_table_ptrs.push_back(etable);
	}

	/*
	* r_vals: the fixed-length Ri identifiers for the combination of each block i
	* i_vals: the variable-length Ii identifiers for the permutation of each block i
	* l_partial_sums: the sum of I encoding lengths up to block i
	* n_partial_sums: the sum of character counts for each character up to block i
	*/
	vector<unsigned int> r_vals; //Note: using unsigned int wastes space!
	vector<unsigned int> i_vals; //Note: using unsigned int wastes space!
	vector<unsigned int> l_partial_sums;
	vector<unsigned int> n_partial_sums ((num_blocks+1)*r, 0); // 2-D: blocks and chars

	l_partial_sums.push_back(0);

	string current_block;
	unsigned int curr_r_index, curr_i_index, curr_l;
	unsigned int block_count = 0;
	for (i=0; i<=n-u; i += u) { //Iterate over every full block
		if (i+u > n) current_block = s.substr(i);
		else current_block = s.substr(i, u);

		// Inefficient: we don't really care about compression efficiency, but this could be better
		// Figure out which index in R this is
		vector<unsigned int> curr_combination (r,0); // Tranform string into character vector
		for (j=0; j<current_block.length(); j++) {
			for (k=0; k<r; k++) {
				if (current_block.at(j) == nseq_alphabet.at(k)) {
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
				r_vals.push_back(j);
				curr_r_index = j;
				break;
			}
		}

		// Next, find permutation Ii
		E_table* table = E_table_ptrs.at(curr_r_index);
		curr_l = ceil(log(table->entries.size()));
		l_partial_sums.push_back(l_partial_sums.back() + curr_l);
		for (j=0; j<table->entries.size(); j++) {
			G_entry* entry = table->entries.at(j);
			if(current_block.compare(entry->sequence) == 0) {
				i_vals.push_back(j);
				// update n_partial_sums:
				for (k=0; k<r; k++) {
					n_partial_sums.at((block_count+1)*r+k) = n_partial_sums.at(block_count*r+k) + curr_combination.at(k);
				}
				break;
			}
		}
		block_count++;
	}
	string rmdr = s.substr(i);
/*
	// DEBUG
	printf("n = %i; r = %i; u = %i\n", (int)n, (int)r, (int)u);
	printf ("Vector size is %u\n", (unsigned int) combinations.size());
	printf("Num blocks = %u\n", (unsigned int) num_blocks);
	unsigned int count = 0;
	unsigned int myindex = 0;
	for(vector<vector<unsigned int> >::iterator it = combinations.begin(); it != combinations.end(); ++it) {
		printf ("%u: (", myindex);
		for (count=0; count<it->size(); count++) {
			printf("%u",it->at(count));
			if (count != it->size() - 1) printf(",");
		}
		printf(")\n");
		myindex++;
	}
	for (i=0; i<combinations.size(); i++) {
		E_table* etable = E_table_ptrs.at(i);
		
		printf("\nE Table for i = %u:\n", (unsigned int) i);
		printf("------------------\n");
		for (unsigned int j=0; j<etable->entries.size(); j++) {
			G_entry* entry = etable->entries.at(j);
			printf("%u: %s\n", j, entry->sequence.c_str());
			for(unsigned int k=0; k<r; k++) {
				for(unsigned int l=0; l<u; l++) {
					printf("  %i", entry->ranks.at(k*u+l));
				}
				printf("\n");
			}
		}
	}
	printf("R Vals:\n-----------\n");
	for (i=0; i<r_vals.size(); i++) {
		printf("%u (", r_vals.at(i));
		for (j=0; j<r; j++) {
			printf("%u", combinations.at(r_vals.at(i)).at(j));
			if (j != r-1) printf(",");
		}
		printf(")\n");
	}
	printf("L partial sums:\n-------------\n");
	for (i=0; i<l_partial_sums.size(); i++) {
		printf("%u ", l_partial_sums.at(i));
	}
	printf("\n");
	printf("N partial sums:\n-------------\n");
	for (i=0; i<num_blocks; i++) {
		for (j=0; j<r; j++) {
			printf("%u ", n_partial_sums.at(i*r+j));
		}
		printf("\n");
	}
	printf("\n");
	// END DEBUGGING
*/
	printf("Bit sequences:\n");
	
	printf("E table:\n-----------\n");
		for (i=0; i<combinations.size(); i++) {
		E_table* etable = E_table_ptrs.at(i);
		
		for (unsigned int j=0; j<etable->entries.size(); j++) {
			G_entry* entry = etable->entries.at(j);
			cout << j;
			for(unsigned int k=0; k<r; k++) {
				for(unsigned int l=0; l<u; l++) {
				  cout << entry->ranks.at(k*u+l);
				}
				printf("\n");
			}
		}
	}
		printf("\n");
	printf("I Vals:\n-----------\n");
	for (i=0; i<i_vals.size(); i++) {
	  cout << toBinStr(i_vals.at(i));
		for (j=0; j<r; j++) {
		  cout << toBinStr(combinations.at(i_vals.at(i)).at(j));
		}
	}
	printf("\n");
	printf("R Vals:\n-----------\n");
	for (i=0; i<r_vals.size(); i++) {
	  cout << toBinStr(r_vals.at(i));
		for (j=0; j<r; j++) {
		  cout << toBinStr(combinations.at(r_vals.at(i)).at(j));
		}
	}
	printf("\n");
	printf("L partial sums:\n-------------\n");
	for (i=0; i<l_partial_sums.size(); i++) {
	  cout << toBinStr(l_partial_sums.at(i));
	}

	printf("\n");
	printf("N partial sums:\n-------------\n");
	for (i=0; i<num_blocks; i++) {
		for (j=0; j<r; j++) {
		  cout << toBinStr(n_partial_sums.at(i*r+j));
		}
	}
	printf("\n");

	decompress(E_table_ptrs, r_vals, i_vals, l_partial_sums, n_partial_sums, rmdr);
	char acc = access(u, E_table_ptrs, r_vals, i_vals, l_partial_sums, n_partial_sums,rmdr, 31);
	cout << "31st index = " << acc << endl;

	return "Hello World!"; //Replace "Hello World" with compressed string here
}

string NavarroSeq::decompress(vector<E_table*>& E_table_ptrs, vector<unsigned int>& r_vals, vector<unsigned int>& i_vals, vector<unsigned int>& l_partial_sums, vector<unsigned int>& n_partial_sums, string rmdr)
{
	// For now, just assume we already have these to work with.
	//vector<E_table*> E_table_ptrs;
	//vector<unsigned int> r_vals; //Note: using unsigned int wastes space!
	//vector<unsigned int> i_vals; //Note: using unsigned int wastes space!
	//vector<unsigned int> l_partial_sums;
	//vector<unsigned int> n_partial_sums; // 2-D: blocks and chars

	string s = "";

	for (unsigned int i=0; i<r_vals.size(); i++) {
		E_table* table = E_table_ptrs.at(r_vals.at(i));
		G_entry* entry = table->entries.at(i_vals.at(i));
		s += entry->sequence;
	}

	cout << "DECOMPRESSED STRING = " << s << rmdr << endl;
	return "Hello again, world!";
}

char NavarroSeq::access(unsigned int u, vector<E_table*>& E_table_ptrs, vector<unsigned int>& r_vals, vector<unsigned int>& i_vals, vector<unsigned int>& l_partial_sums, vector<unsigned int>& n_partial_sums, string rmdr, int index)
{
	unsigned int block = floor(index/u);
	unsigned int l = index - block*u;
	E_table* table = E_table_ptrs.at(r_vals.at(block));
	G_entry* entry = table->entries.at(i_vals.at(block));
	char retval = entry->sequence.at(l);
	return retval;
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
			vector<unsigned int> values_copy (values);
			values_copy.at(i)++;
			enumerate(values_copy, combos, i, r, u-1);
		}
	}
}

void NavarroSeq::get_etable_rows(string prefix, vector<char> ranks, vector<unsigned int> combination, E_table* table)
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
		string next_prefix = prefix + nseq_alphabet.at(i);
		vector<unsigned int> remaining_combo (combination);
		remaining_combo.at(i)--;
		for (j=i*nseq_u + prefix.length(); j<(i+1)*nseq_u; j++) {
			ranks.at(j)++;
		}
		get_etable_rows(next_prefix, ranks, remaining_combo, table);
	}
}

E_table* NavarroSeq::get_etable(vector<unsigned int> combination)
{
	E_table* table = new E_table();
	unsigned int i,j;
	for (i=0; i<combination.size(); i++) {
		if(combination.at(i) == 0) continue;
		string prefix = "";
		prefix += nseq_alphabet.at(i);
		vector<char> ranks (nseq_r*nseq_u, 0);
		vector<unsigned int> remaining_combo (combination);
		remaining_combo.at(i)--;
		for (j=i*nseq_u; j<(i+1)*nseq_u; j++) {
			ranks.at(j)++;
		}
		get_etable_rows(prefix, ranks, remaining_combo, table);
	}
	return table;
}

int main() {
	// test sequence
	NavarroSeq::compress("abbabaababbababcabbabaababbababcabbabaababbababcabbabaababbababcabbabaababbababcabc");

	return 0;
}
