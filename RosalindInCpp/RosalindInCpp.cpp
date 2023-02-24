#include <iostream>
#include <string>
#include <fstream>
#include <unordered_map>
#include "FastaTools.h"

using namespace std;

int main()
{
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

    return 0;
}

