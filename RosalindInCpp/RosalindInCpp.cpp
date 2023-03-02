#include <iostream>
#include <string>
#include <fstream>
#include <unordered_map>
#include "FastaTools.h"

using namespace std;

void TestTool();
unordered_map<string, string> ReadCodonTable();

//read codons - will be accessible by any instances of FastaTools
unordered_map<string, string> codon_table = ReadCodonTable();

int main()
{
    TestTool();
    return 0;
}

// test tools by calling them and comparing against the Python code which is already verified.
void TestTool() {

    string f_path;

    cout << "File path (without \"\"): " << endl; //ask for file path
    getline(cin, f_path);
    ifstream checkFile(f_path);

    while (!checkFile.is_open()) { //test if file path is correct
        cout << "Error with file path: " << f_path << "\nEnter new file path:\n";
        getline(cin, f_path);
        checkFile.open(f_path);
    }
    checkFile.close();

    //Create an instance of ReadFasta class and print the sequence ids and the sequences
    FastaTools fasta(f_path);
    vector<string> seqs = fasta.GetFastaSeqs();
    string seq = seqs[0];
    unordered_map<char, int> counts = fasta.countNucs(seq);

    for (auto i : counts)
        cout << i.first << ": " << i.second << endl;

    cout << "Original Sequence:\t\t5' " << seq << " 3'" << endl;

    cout << "Transcribed Sequence:\t\t5' " << fasta.TranscribeDNA(seq) << " 3'" << endl;

    cout << "Reverse complement Sequence:\t3' " << fasta.RevComplDNA(seq) << " 5'" << endl;

    pair<string, float> gc_content = fasta.HighestGC();
    cout << "Highest GC content in sequence with ID: " << gc_content.first << endl << "with a content of: " << gc_content.second << "%" << endl;

    cout << "Hamming Distance: " << fasta.HammingDist(seqs[0], seqs[1]) << endl;

    cout << "Translated Protein: " << fasta.TranslateRNA(fasta.TranscribeDNA(seq), codon_table) << endl;

    auto find_motifs = fasta.FindMotif("GTGCTCTCTCTTAACTCTTAACGTCTTAACTCTTAACTCTTAACCTATTCTTAACTATCTTAACTTGTCTTAACCTCTTAACTCTTAACTCTTAACCTCTTAACTGGAGCGTCTCTTAACACGGCGTCTTAACGGATTCTTAACAAACTATCTTAACTCTTAACCTCTTAACCGTTCTTAACATCTTAACTTTCTTAACAATCTTAACTCTTAACCTCTTAACAATCTTAACTCTCTTAACGTCTTAACGAATCTTAACGATATCTTAACATCTTTCTATAGAATAGTGTCTTAACTCCCGGCTCTTAACTTCTTAACGTCTTAACGTGTCTTAACTTTCTTAACCTCTTAACGGGTCTTAACTACTACTCTTAACTCTTAACAGTTACCGTCTTAACTCTCTTAACACTTCTTAACTCTTAACTTCTTAACTTTCTTAACTCTTAACTCTTAACTAATCCTCTTAACATTCTTAACTTACTTCGTTTCTTAACTCTTAACTCCATCTTAACAATCTTAACTGCTCTTAACGTCTTAACACTCTTAACTCTTAACGGTGAGATCTTAACTGATCTTAACGTCTTAACATCTTAACTGCTCAGTCTTAACTCTTAACTCTTAACGTCTTAACCTCTTAACGGTCTTAACGTGACAATCAGTCTTAACTCTTAACTTTCTTAACTGTTCTTAACATCTTAACATTCTTAACTCTCCTCTTAACGAGTCTTAACTCTTAACTCTTAACCATTTTACCTCCGGATAGGTACGATCTTAACATCTTAACTCTTAACTCTTAACTCTTAACTGTCTTAACTCTTAACTTCTTAACACTCTTAACTTCTTCCTTCTTAACTCTTAACTGTTCTTAAC", "TCTTAACTC");
    cout << "Motif found at position(s):\n"; 
    for (int i = 0; i < find_motifs.size(); i++)
        cout << find_motifs[i] << " ";
}

unordered_map<string, string> ReadCodonTable() {
    unordered_map<string, string> codons;
    string current_line;
    ifstream readFileStream("../rosalind_exercises_python/codon_tab.txt");

    while (getline(readFileStream, current_line)) {
        codons[current_line.substr(0, 3)] = current_line.substr(4, 4);
    }
    readFileStream.close();

    return codons;
}