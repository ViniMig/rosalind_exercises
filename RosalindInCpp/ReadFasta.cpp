#include <iostream>
#include <string>
#include <string.h>
#include <unordered_map>
#include <fstream>
#include "ReadFasta.h"

using namespace std;

//constructor
ReadFasta::ReadFasta(string fpath)
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
	readFileStream.close();

	fastaHMap = fastaMap;
}

//methods
void ReadFasta::ShowFastaIds() {
	for (auto i = fastaHMap.begin(); i != fastaHMap.end(); i++) {
		cout << i->first << ";\n";
	}
}

void ReadFasta::ShowFastaSeqs() {
	for (auto i = fastaHMap.begin(); i != fastaHMap.end(); i++) {
		cout << i->second << ";\n";
	}
}

vector<string> ReadFasta::GetFastaIds() {

	vector<string> fastaIds(fastaHMap.size());

	for (auto i = fastaHMap.begin(); i != fastaHMap.end(); i++) {
		fastaIds.push_back(i->first);
	}

	return fastaIds;
}

vector<string> ReadFasta::GetFastaSeqs() {

	vector<string> fastaSeqs(fastaHMap.size());

	for (auto i = fastaHMap.begin(); i != fastaHMap.end(); i++) {
		fastaSeqs.push_back(i->second);
	}

	return fastaSeqs;
}