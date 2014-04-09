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

vector<vector<unsigned int> > get_frequencies(string filename)
{
  ifstream file(filename);
  if(!file)
    {
      cout << "Could not find file." << endl;
    }

  vector<vector<unsigned int> > frequencies;
  
for(int i = 0; i < 9; i++)
    {
      //First 3 correspond to zero, second 3 correspond to 1, third 3 correspond to 2
      vector<unsigned int> values;
      frequencies.push_back(values);
    }

  char last = 'A';
  char current;
  string window;

  for(unsigned int i = 0; i < get_size(filename); i++)
    {
      current = file.get();
      //for test data
      if(current == '\n' || current == ' ')
	{
	  i++;
	  current = file.get();
	}

      if(i == get_size(filename))
	break;

      //cout << current;
      if(current == '0' && (last == '1' || last == '2'))
	{
	  if(last == '1')
	    {
	      while(frequencies.at(3).size() <= window.size()+1)
	       {
		 frequencies.at(3).push_back(0);
	       }
	      frequencies.at(3).at(window.size())++;
	    } 
	  else if(last == '2')
	    {
	      while(frequencies.at(6).size() <= window.size()+1)
		{
		  frequencies.at(6).push_back(0);
		}	 
	      frequencies.at(6).at(window.size())++;
	    }
	  window.clear();
	}
      else if((current == '1' || current == '2') && last == '0')
	{
	  if(current == '1')
	    {
	    while(frequencies.at(1).size() <= window.size()+1)
		{
		  frequencies.at(1).push_back(0);
		}	 
	      frequencies.at(1).at(window.size())++;
	    }
	  else if(current == '2')
	    {
	      while(frequencies.at(2).size() <= window.size()+1)
		{
		  frequencies.at(2).push_back(0);
		}	 
	      frequencies.at(2).at(window.size())++;
	    }
	  window.clear();
	}

      //add current element to window
      window.append(1,current);

      //check to see if it's 0*
      if(regex_match(window, zero))
	{
	  if(current == '0')
	    {
	      while(frequencies.at(0).size() <= window.size())
		{
		  frequencies.at(0).push_back(0);
		}	 
	      frequencies.at(0).at(window.size()-1)++;
	    }
	  else
	    { cout << "error!" << endl; }
	}

      //check to see if it's {1,2}*1
      else if(regex_match(window, one))
	{
	  if(current == '1')
	    {
	      while(frequencies.at(4).size() <= window.size())
		{
		  frequencies.at(4).push_back(0);
		}	 
	      frequencies.at(4).at(window.size()-1)++;
	    }
	  else if(current == '2')
	    {
	      //actually unused
	      while(frequencies.at(5).size() <= window.size())
		{
		  frequencies.at(5).push_back(0);
		}	 
	      frequencies.at(5).at(window.size()-1)++;
	    }
	}

      //check to see if it's {1,2}*2
      else if(regex_match(window, two))
	{
	  if(current == '1')
	    {
	      //actually unused
	      while(frequencies.at(7).size() <= window.size())
		{
		  frequencies.at(7).push_back(0);
		}	 
	      frequencies.at(7).at(window.size()-1)++;
	    }
	  else if(current == '2')
	    {
	      while(frequencies.at(8).size() <= window.size())
		{
		  frequencies.at(8).push_back(0);
		}	 
	      frequencies.at(8).at(window.size()-1)++;
	    }
	}
      else
	cout << "error: " << window << endl;
      last = current;
    }
  //cout << endl;
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
      else if(regex_match(window, one) || regex_match(window, two))
	{
	  int depth = window.size();
	  int temp1Zero, temp2Zero, tempOne, tempTwo;
	  if(depth < frequencies.at(3).size())
	    temp1Zero = frequencies.at(3).at(depth);
	  else
	    temp1Zero = frequencies.at(3).at(frequencies.at(3).size()-1);

	  if(depth < frequencies.at(4).size())
	    tempOne = frequencies.at(4).at(depth);
	  else
	    tempOne = frequencies.at(4).at(frequencies.at(4).size()-1);

	  if(depth < frequencies.at(6).size())
	    temp2Zero = frequencies.at(6).at(depth);
	  else
	    temp2Zero = frequencies.at(6).at(frequencies.at(6).size()-1);

	  if(depth < frequencies.at(8).size())
	    tempTwo = frequencies.at(8).at(depth);
	  else
	    tempTwo = frequencies.at(8).at(frequencies.at(8).size()-1);

	  int tempZero = tempOne + tempTwo;
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
  if(argc < 2)
    {
      cout << "Usage: " << argv[0] << "inputfilename sizefile numbertestfiles" << endl;
      return 0;
    }

  string filename = argv[1];
  unsigned int size = atoi(argv[2]);
  unsigned int fileNum = atoi(argv[3]);

  /* Playing around with regex
  string s = "11112";
  if(regex_match(s,zero))
    cout << "0" << endl; 
  else if(regex_match(s,one))
	  cout << "1" << endl;
  else if(regex_match(s,two))
  cout << "2" << endl; */
  cout << "Getting info..." << endl;
  vector<vector<unsigned int> > freq =  get_frequencies(filename);
  cout << "Generating test data..." << endl;
  generate_data("outtest.txt", size, freq);

  /* Testing
  ifstream ifs("outtest.txt");
  while(ifs.good())
    {
      char x = ifs.get();
      if(ifs.good())
	cout << x;
    }
  */
}
