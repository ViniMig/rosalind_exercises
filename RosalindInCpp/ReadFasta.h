#include <iostream>
#include <string>
#include <unordered_map>

using namespace std;

class ReadFasta {
	private:
		unordered_map<string, string> fastaHMap;

	public:
		ReadFasta(string fpath);
		void ShowFastaIds();
		void ShowFastaSeqs();
		vector<string> GetFastaIds();
		vector<string> GetFastaSeqs();
};