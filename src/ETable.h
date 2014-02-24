#ifndef E_table
#define ETABLE

using namespace std;

struct G_entry {
	string sequence;
	vector<char> ranks;

	G_entry(string seq, vector<char> ranks_in)
	{
		sequence = seq;
		ranks = ranks_in;
	}
};

struct E_table {
	
	vector<G_entry*> entries;

	void add_entry(G_entry* entry) {
		entries.push_back(entry);
	}
};

#endif