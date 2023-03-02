#include <iostream>
#include <string>
#include <string.h>
#include <unordered_map>
#include <fstream>
#include <regex>
#include "fasta_tools.h"

using namespace std;

//constructor
FastaTools::FastaTools(string f_path)
	:fasta_h_map_{}
{
	unordered_map<string, string> fasta_map;
	int i = 0;
	string current_line;
	string current_seq = "";
	string current_id;
	bool is_start = true;

	ifstream readFileStream(f_path);

	while (getline(readFileStream, current_line)) {
		if (current_line[0] == '>') {
			if (is_start) {
				is_start = false;
			}
			else {
				fasta_map[current_id] = current_seq;
				current_seq = "";
			}
			i = 0;
			current_id = current_line.substr(1, current_line.size() - 1);
		}
		else {
			current_seq += current_line;
			i++;
		}
	}
	fasta_map[current_id] = current_seq; //add the last id and sequence of the file
	readFileStream.close();

	fasta_h_map_ = fasta_map;
}

//methods

// Display IDs of the fasta sequences on the console
void FastaTools::ShowFastaIds() {
	for (auto i = fasta_h_map_.begin(); i != fasta_h_map_.end(); i++) {
		cout << i->first << ";\n";
	}
}

// Display sequences present in the fasta file on the console
void FastaTools::ShowFastaSeqs() {
	for (auto i = fasta_h_map_.begin(); i != fasta_h_map_.end(); i++) {
		cout << i->second << ";\n";
	}
}

// Return a vector of the IDs
vector<string> FastaTools::GetFastaIds() {

	vector<string> fasta_ids;

	for (auto i = fasta_h_map_.begin(); i != fasta_h_map_.end(); i++) {
		fasta_ids.push_back(i->first);
	}

	return fasta_ids;
}

// Return a vector of the Sequences - can be used for further analysis.
vector<string> FastaTools::GetFastaSeqs() {

	vector<string> fasta_seqs;

	for (auto i = fasta_h_map_.begin(); i != fasta_h_map_.end(); i++) {
		fasta_seqs.push_back(i->second);
	}

	return fasta_seqs;
}

// Return a dictionary with the counts as values for each of the 4 nucleotides present.
unordered_map<char, int> FastaTools::countNucs(string seq) {
	unordered_map<char, int> nuc_count;
	nuc_count['A'] = 0;
	nuc_count['C'] = 0;
	nuc_count['T'] = 0;
	nuc_count['G'] = 0;

	for (char nuc : seq) {
		if (nuc == 'A')
			nuc_count['A']++;
		else if (nuc == 'C')
			nuc_count['C']++;
		else if (nuc == 'T')
			nuc_count['T']++;
		else if (nuc == 'G')
			nuc_count['G']++;
		else {
			cout << "Invalid nucleotide. This is an invalid sequence, check for errors!" << endl;
			break;
		}
	}

	return nuc_count;
}

// Return the RNA sequence string of the transcribed DNA.
string FastaTools::TranscribeDNA(string seq) {
	string transcr_seq = "";

	for (char nuc : seq) {
		if (nuc == 'T')
			nuc = 'U';

		transcr_seq += nuc;
	}

	return transcr_seq;
}

// Return the reverse complement sequence of the input sequence.
string FastaTools::RevComplDNA(string seq) {
	string rev_compl_seq = "";

	for (int i = seq.size() - 1; i >= 0; i--) {
		switch (seq[i]) {
		case 'A':
			rev_compl_seq += 'T';
			break;
		case 'T':
			rev_compl_seq += 'A';
			break;
		case 'C':
			rev_compl_seq += 'G';
			break;
		case 'G':
			rev_compl_seq += 'C';
			break;
		default:
			return "Invalid Sequence! Non nucleotide present: " + seq[i];
		}
	}

	return rev_compl_seq;
}

//Return ID of the sequence with the highest GC content in te fasta file, and the content in %.
pair<string, float> FastaTools::HighestGC() {
	pair <string, float> high_GC_cont;
	high_GC_cont.first = "";
	high_GC_cont.second = 0;
	vector<string> fasta_ids = GetFastaIds();
	vector<string> sequences = GetFastaSeqs();

	for (int i = 0; i < sequences.size(); i++) {
		unordered_map<char, int> nuc_count = countNucs(sequences[i]);
		float curr_content = ((nuc_count['C'] + nuc_count['G']) / (sequences[i].size() * 1.0F)) * 100;
		if (curr_content > high_GC_cont.second) {
			high_GC_cont.first = fasta_ids[i];
			high_GC_cont.second = curr_content;
		}
	}
	return high_GC_cont;
}

// Return the hamming distance between any 2 sequences of same size.
int FastaTools::HammingDist(string seq_1, string seq_2) {
	int hamm_dist = 0;

	if (seq_1.size() != seq_2.size()) {
		cout << "Error! Sequences with different sizes!" << endl;
		return -1;
	}

	for (int i = 0; i < seq_1.size(); i++) {
		if (seq_1[i] != seq_2[i])
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

//Return list of positions where motif is present in the sequence
vector<int> FastaTools::FindMotif(string seq, string motif) {
	vector<int> motif_pos;

	string::size_type start_pos = 0;
	while (string::npos != (start_pos = seq.find(motif, start_pos)))
	{
		motif_pos.push_back(start_pos + 1); //adding one because of 0 index.
		start_pos++;
	}

	return motif_pos;
}