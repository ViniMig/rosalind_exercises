#include <iostream>
#include <string>
#include <fstream>
#include <unordered_map>
#include "FastaTools.h"

using namespace std;

void TestTool();

int main()
{
    TestTool();
    return 0;
}

// test tools by calling them and comparing against the Python code which is already verified.
void TestTool() {
    string fPath;

    cout << "File path (without \"\"): " << endl; //ask for file path
    getline(cin, fPath);
    ifstream checkFile(fPath);

    while (!checkFile.is_open()) { //test if file path is correct
        cout << "Error with file path: " << fPath << "\nEnter new file path:\n";
        getline(cin, fPath);
        checkFile.open(fPath);
    }
    checkFile.close();

    //Create an instance of ReadFasta class and print the sequence ids and the sequences
    FastaTools fasta(fPath);
    vector<string> seqs = fasta.GetFastaSeqs();
    string seq = seqs[0];
    unordered_map<char, int> counts = fasta.countNucs(seq);

    for (auto i : counts)
        cout << i.first << ": " << i.second << endl;

    cout << "Original Sequence:\t\t5' " << seq << " 3'" << endl;

    cout << "Transcribed Sequence:\t\t5' " << fasta.transcribeDNA(seq) << " 3'" << endl;

    cout << "Reverse complement Sequence:\t3' " << fasta.revComplDNA(seq) << " 5'" << endl;

    pair<string, float> gc_content = fasta.HighestGC();
    cout << "Highest GC content in sequence with ID: " << gc_content.first << endl << "with a content of: " << gc_content.second << "%" << endl;

    cout << "Hamming Distance: " << fasta.HammingDist(seqs[0], seqs[1]) << endl;
}