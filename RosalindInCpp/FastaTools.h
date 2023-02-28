#include <iostream>
#include <string>
#include <unordered_map>

using namespace std;

class FastaTools {
	private:
		unordered_map<string, string> fastaHMap;

	public:
		FastaTools(string fpath);
		void ShowFastaIds();
		void ShowFastaSeqs();
		vector<string> GetFastaIds();
		vector<string> GetFastaSeqs();
		unordered_map<char, int> countNucs(string seq);
		string TranscribeDNA(string seq);
		string RevComplDNA(string seq);
		pair<string, float> HighestGC();
		int HammingDist(string seq1, string seq2);
		string TranslateRNA(string seq, unordered_map<string, string> codons);
		vector<int> FindMotif(string seq, string motif);
};