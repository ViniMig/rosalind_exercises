#include <iostream>
#include <string>
#include <string.h>
#include <unordered_map>
#include <fstream>
#include "FastaTools.h"

using namespace std;

//constructor
FastaTools::FastaTools(string fpath)
	:fastaHMap()
{
	unordered_map<string, string> fastaMap;
	int i = 0;
	string currentLine;
	string currentSeq = "";
	string currentID;
	bool isStart = true;

	ifstream readFileStream(fpath);

	while (getline(readFileStream, currentLine)) {
		if (currentLine[0] == '>') {
			if (isStart) {
				isStart = false;
			}
			else {
				fastaMap[currentID] = currentSeq;
				currentSeq = "";
			}
			i = 0;
			currentID = currentLine.substr(1, currentLine.size() - 1);
		}
		else {
			currentSeq += currentLine;
			i++;
		}
	}
	fastaMap[currentID] = currentSeq; //add the last id and sequence of the file
	readFileStream.close();

	fastaHMap = fastaMap;
}

//methods
void FastaTools::ShowFastaIds() {
	for (auto i = fastaHMap.begin(); i != fastaHMap.end(); i++) {
		cout << i->first << ";\n";
	}
}

void FastaTools::ShowFastaSeqs() {
	for (auto i = fastaHMap.begin(); i != fastaHMap.end(); i++) {
		cout << i->second << ";\n";
	}
}

vector<string> FastaTools::GetFastaIds() {

	vector<string> fastaIds(fastaHMap.size());

	for (auto i = fastaHMap.begin(); i != fastaHMap.end(); i++) {
		fastaIds.push_back(i->first);
	}

	return fastaIds;
}

vector<string> FastaTools::GetFastaSeqs() {

	vector<string> fastaSeqs;

	for (auto i = fastaHMap.begin(); i != fastaHMap.end(); i++) {
		fastaSeqs.push_back(i->second);
	}

	return fastaSeqs;
}

unordered_map<char, int> FastaTools::countNucs(string seq) {
	unordered_map<char, int> nucCount;
	nucCount['A'] = 0;
	nucCount['C'] = 0;
	nucCount['T'] = 0;
	nucCount['G'] = 0;

	cout	<< "Sequence: " << seq << endl
			<< "Freqs:\n";
	for (char nuc : seq) {
		if (nuc == 'A')
			nucCount['A']++;
		else if (nuc == 'C')
			nucCount['C']++;
		else if (nuc == 'T')
			nucCount['T']++;
		else if (nuc == 'G')
			nucCount['G']++;
		else {
			cout << "Invalid nucleotide. This is an invalid sequence, check for errors!" << endl;
			break;
		}
	}

	return nucCount;
}