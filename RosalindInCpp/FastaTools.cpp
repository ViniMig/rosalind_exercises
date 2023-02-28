#include <iostream>
#include <string>
#include <string.h>
#include <unordered_map>
#include <fstream>
#include "FastaTools.h"

using namespace std;

//constructor
FastaTools::FastaTools(string fpath)
	:fastaHMap{}
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

// Display IDs of the fasta sequences on the console
void FastaTools::ShowFastaIds() {
	for (auto i = fastaHMap.begin(); i != fastaHMap.end(); i++) {
		cout << i->first << ";\n";
	}
}

// Display sequences present in the fasta file on the console
void FastaTools::ShowFastaSeqs() {
	for (auto i = fastaHMap.begin(); i != fastaHMap.end(); i++) {
		cout << i->second << ";\n";
	}
}

// Return a vector of the IDs
vector<string> FastaTools::GetFastaIds() {

	vector<string> fastaIds;

	for (auto i = fastaHMap.begin(); i != fastaHMap.end(); i++) {
		fastaIds.push_back(i->first);
	}

	return fastaIds;
}

// Return a vector of the Sequences - can be used for further analysis.
vector<string> FastaTools::GetFastaSeqs() {

	vector<string> fastaSeqs;

	for (auto i = fastaHMap.begin(); i != fastaHMap.end(); i++) {
		fastaSeqs.push_back(i->second);
	}

	return fastaSeqs;
}

// Return a dictionary with the counts as values for each of the 4 nucleotides present.
unordered_map<char, int> FastaTools::countNucs(string seq) {
	unordered_map<char, int> nucCount;
	nucCount['A'] = 0;
	nucCount['C'] = 0;
	nucCount['T'] = 0;
	nucCount['G'] = 0;

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

// Return the RNA sequence string of the transcribed DNA.
string FastaTools::TranscribeDNA(string seq) {
	string tSeq = "";

	for (char nuc : seq) {
		if (nuc == 'T')
			nuc = 'U';

		tSeq += nuc;
	}

	return tSeq;
}

// Return the reverse complement sequence of the input sequence.
string FastaTools::RevComplDNA(string seq) {
	string revCompSeq = "";

	for (int i = seq.size() - 1; i >= 0; i--) {
		switch (seq[i]) {
		case 'A':
			revCompSeq += 'T';
			break;
		case 'T':
			revCompSeq += 'A';
			break;
		case 'C':
			revCompSeq += 'G';
			break;
		case 'G':
			revCompSeq += 'C';
			break;
		default:
			return "Invalid Sequence! Non nucleotide present: " + seq[i];
		}
	}

	return revCompSeq;
}

//Return ID of the sequence with the highest GC content in te fasta file, and the content in %.
pair<string, float> FastaTools::HighestGC() {
	pair <string, float> highestGC;
	highestGC.first = "";
	highestGC.second = 0;
	vector<string> fastaIds = GetFastaIds();
	vector<string> sequences = GetFastaSeqs();

	for (int i = 0; i < sequences.size(); i++) {
		unordered_map<char, int> nucCount = countNucs(sequences[i]);
		float currContent = ((nucCount['C'] + nucCount['G']) / (sequences[i].size() * 1.0F)) * 100;
		if (currContent > highestGC.second) {
			highestGC.first = fastaIds[i];
			highestGC.second = currContent;
		}
	}
	return highestGC;
}

// Return the hamming distance between any 2 sequences of same size.
int FastaTools::HammingDist(string seq1, string seq2) {
	int hamm_dist = 0;

	if (seq1.size() != seq2.size()) {
		cout << "Error! Sequences with different sizes!" << endl;
		return -1;
	}

	for (int i = 0; i < seq1.size(); i++) {
		if (seq1[i] != seq2[i])
			hamm_dist++;
	}

	return hamm_dist;
}

// Return the translated sequence until the stop codon in a given sequence.
string FastaTools::TranslateRNA(string seq, unordered_map<string, string> codons) {
	string rna_seq = "";

	for (int i = 0; i < seq.size(); i+=3) {
		string amino_acid = codons[seq.substr(i, 3)];
		if (amino_acid == "Stop")
			break;

		rna_seq += amino_acid;
	}
	
	return rna_seq;
}