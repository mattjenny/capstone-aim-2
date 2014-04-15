#include <regex>
#include <string.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <random>

using namespace std;
regex zero("(0*)");
regex one("(1|2)*(1)");
regex two("(1|2)*(2)");

unsigned int get_size(string filename)
{
  ifstream file(filename);
  file.seekg(0, file.end);
  int size = file.tellg();
  return size;
}

vector<vector<unsigned int> > get_freq_file(string filename)
{
  ifstream file(filename);
  if(!file)
    cout << "Could not find " << filename << endl;

  vector<vector<unsigned int> > frequencies;

  for(int i = 0; i < 9; i++)
    {
      //First 3 correspond to zero, second 3 correspond to 1, third 3 correspond to 2
      vector<unsigned int> values(20, 0);
      frequencies.push_back(values);
    }

  string line;
  int i = 0;
  while(getline(file,line))
    {
      stringstream ss(line);

      int j = 0;
      string x, throwaway;
      getline(ss, throwaway, ' ');
      while(ss)
	{
	  if(j > 19)
	    break;
	  getline(ss, x, ' ');
	  frequencies.at(i).at(j) = stoi(x);
	  j++;
	}
      i++;
    }

  return frequencies;
}

vector<vector<unsigned int> > get_frequencies(string filename)
{
  ifstream file(filename);
  if(!file)
      cout << "Could not find " << filename << endl;

  string frequencies_out = "frequencies.txt";
  ofstream ofs(frequencies_out);

  vector<vector<unsigned int> > frequencies;
  
for(int i = 0; i < 9; i++)
    {
      //First 3 correspond to zero, second 3 correspond to 1, third 3 correspond to 2
      vector<unsigned int> values(20, 0);
      frequencies.push_back(values);
    }

  char last = 'A';
  char current;
  string window;
  cout << "file size: " << get_size(filename) << endl;

  int max = get_size(filename);
  //set bounds on filesize if needed
  /*int max = 1000000;
  if(get_size(filename)<max)
  max = get_size(filename);*/
  int div = max/100;
  int count = div;

  for(unsigned int i = 0; i < max; i++)
    {
      current = file.get();
      //for test data
      if(current == '\n' || current == ' ')
	{
	  i++;
	  current = file.get();
	}

      if(i == max)
	break;

      //cout << current;
      
      if(i == count)
	{
	  cout << count << endl;
	  count = count+div;
	}
      //check to see if it's 0*
      if(regex_match(window, zero))
	{
	  if(current == '0')
	    {
	      int index = window.size();
	      if(index>19)
		index = 19;
	      frequencies.at(0).at(index)++;
	    }
	  else if(current == '1')
	    {
	      int index = window.size();
	      if(index>19)
		index = 19;	 
	      frequencies.at(1).at(index)++;
	    }
	  else if(current == '2')
	    {
	      int index = window.size();
	      if(index>19)
		index = 19;	 
	      frequencies.at(2).at(index)++;
	    }
	}

      //check to see if it's {1,2}*1
      else if(regex_match(window, one))
	{
	  if(current == '0')
	    {
	      int index = window.size();
	      if(index>19)
		index = 19;
	      frequencies.at(3).at(index)++;
	    } 
	  else if(current == '1')
	    {
	      int index = window.size();
	      if(index>19)
		index = 19;	 
	      frequencies.at(4).at(index)++;
	    }
	  else if(current == '2')
	    {
	      int index = window.size();
	      if(index>19)
		index = 19;	 
	      frequencies.at(5).at(index)++;
	    }
	}

      //check to see if it's {1,2}*2
      else if(regex_match(window, two))
	{
	  if(current == '0')
	    {
	      int index = window.size();
	      if(index>19)
		index = 19;	 
	      frequencies.at(6).at(index)++;
	    }
	  else if(current == '1')
	    {
	      int index = window.size();
	      if(index>19)
		index = 19;	 
	      frequencies.at(7).at(index)++;
	    }
	  else if(current == '2')
	    {
	      int index = window.size();
	      if(index>19)
		index = 19;	 
	      frequencies.at(8).at(index)++;
	    }
	}
      else
	cout << "error: " << window << endl;

      //check for window change
      if(current == '0' && (last == '1' || last == '2'))
	{
	  window.clear();
	}

      if((current == '1' || current == '2') && last == '0')
	{
	  window.clear();
	}

      //add current element to window
      window.append(1,current);

      last = current;
    }
  //cout << endl;

  for(int i = 0; i < frequencies.size(); i++)
    {
      ofs << i << ": ";
      for(int j = 0; j < frequencies.at(i).size(); j++)
	{
	  ofs << frequencies.at(i).at(j) << " ";
	}
      ofs << endl;
    }

  return frequencies;
}

void generate_data(string outfile, unsigned int size, vector<vector<unsigned int> > frequencies)
{
  ofstream ofs(outfile);
  //start with 0 to start, for simplification
  ofs.put(0);

  string window;
  char c;
  window.append(1,'0');

  random_device rdev{};
  default_random_engine generator{rdev()};
  uniform_real_distribution<float> distribution(0.0, 1.0);
  
  for(int i = 0; i < size; i++)
    {
      if(regex_match(window, zero))
	{
	  int depth = window.size();
	  if(depth > 19)
	    depth = 19;
	  int tempZero, tempOne, tempTwo;
	  if(depth < frequencies.at(0).size())
	    tempZero = frequencies.at(0).at(depth);
	  else
	    tempZero = frequencies.at(0).at(frequencies.at(0).size()-1);

	  if(depth < frequencies.at(1).size())
	    tempOne = frequencies.at(1).at(depth);
	  else
	    tempOne = frequencies.at(1).at(frequencies.at(1).size()-1);

	  if(depth < frequencies.at(2).size())
	    tempTwo = frequencies.at(2).at(depth);
	  else
	    tempTwo = frequencies.at(2).at(frequencies.at(2).size()-1);

	  int total = tempZero + tempOne + tempTwo;

	  float rand = distribution(generator)*total;
	  int value = rand;
	  if(value < tempZero)
	    c = '0';
	  if(value < (tempZero+tempOne) && value > tempZero)
	    c = '1';
	  if(value <= total && value > (tempZero+tempOne))
	    c = '2';
	  ofs.put(c);
	}
      else if(regex_match(window, one))
	{
	  int depth = window.size();
	  if(depth > 19)
	    depth = 19;
	  int tempZero, tempOne, tempTwo;
	  if(depth < frequencies.at(3).size())
	    tempZero = frequencies.at(3).at(depth);
	  else
	    tempZero = frequencies.at(3).at(frequencies.at(3).size()-1);

	  if(depth < frequencies.at(4).size())
	    tempOne = frequencies.at(4).at(depth);
	  else
	    tempOne = frequencies.at(4).at(frequencies.at(4).size()-1);

	  if(depth < frequencies.at(5).size())
	    tempTwo = frequencies.at(5).at(depth);
	  else
	    tempTwo = frequencies.at(5).at(frequencies.at(5).size()-1);

	  int total = tempZero + tempOne + tempTwo;

	  float rand = distribution(generator)*total;
	  int value = rand;
	  if(value < tempZero)
	    c = '0';
	  if(value < (tempZero+tempOne) && value > tempZero)
	    c = '1';
	  if(value <= total && value > (tempZero+tempOne))
	    c = '2';
	  ofs.put(c);
	}
      else if(regex_match(window, two))
	{
	  int depth = window.size();
	  if(depth > 19)
	    depth = 19;
	  int tempZero, tempOne, tempTwo;
	  if(depth < frequencies.at(6).size())
	    tempZero = frequencies.at(6).at(depth);
	  else
	    tempZero = frequencies.at(6).at(frequencies.at(6).size()-1);

	  if(depth < frequencies.at(7).size())
	    tempOne = frequencies.at(7).at(depth);
	  else
	    tempOne = frequencies.at(7).at(frequencies.at(7).size()-1);

	  if(depth < frequencies.at(8).size())
	    tempTwo = frequencies.at(8).at(depth);
	  else
	    tempTwo = frequencies.at(8).at(frequencies.at(8).size()-1);

	  int total = tempZero + tempOne + tempTwo;

	  float rand = distribution(generator)*total;
	  int value = rand;
	  if(value < tempZero)
	    c = '0';
	  if(value < (tempZero+tempOne) && value > tempZero)
	    c = '1';
	  if(value <= total && value > (tempZero+tempOne))
	    c = '2';
	  ofs.put(c);
	}

      if(window.size() > 0)
	{
	  if(window.at(window.size()-1) == '0' && c != '0')
	    window.clear();
	  else if(window.at(window.size()-1) != '0' && c == '0')
	    window.clear();
	}
      window.append(1,c);
    }
  /*    
  //Testing purposes
    for(int i = 0; i < frequencies.size(); i++)
    {
      cout << i << ": ";
      for(int j = 0; j < frequencies.at(i).size(); j++)
	{
	  cout << frequencies.at(i).at(j) << " ";
	}
      cout << endl;
    }
  */
}

int main(int argc, char** argv)
{
  if(argc < 4)
    {
      cout << "Usage: " << argv[0] << "inputfilename sizefile numbertestfiles" << endl;
      cout << "Or use -f flag to use frequencies.txt (in place of inputfilename)" << endl;
      return 0;
    }

  bool flag = false;

  string filename = argv[1];
  unsigned int size = atoi(argv[2]);
  unsigned int fileNum = atoi(argv[3]);

  if(strcmp("-f", argv[1]) == 0)
    flag = true;

  cout << "Getting info..." << endl;
  vector<vector<unsigned int> > freq;

  if(flag)
    freq = get_freq_file("frequencies.txt");
  else
    freq = get_frequencies(filename);

  cout << "Generating test data..." << endl;

  string outfileName = "testData";
  string outfileExtension = ".txt";
  for(int i = 1; i <= fileNum; i++)
    {
      stringstream ss;
      ss << i;
      string num = ss.str();
      string outfile = outfileName + num + outfileExtension;
      generate_data(outfile, size, freq);
    }

  /* Testing
  ifstream ifs("outtest.txt");
  while(ifs.good())
    {
      char x = ifs.get();
      if(ifs.good())
	cout << x;
    }
  */

  /* Playing around with regex
  string s = "11112";
  if(regex_match(s,zero))
    cout << "0" << endl; 
  else if(regex_match(s,one))
	  cout << "1" << endl;
  else if(regex_match(s,two))
  cout << "2" << endl; */
}
