/* 
 * File:   main.cc
 * Author: cecco
 *
 * Created on May 29, 2015, 4:18 PM
 */

#include <cstdlib>
#include <SuperImpositor.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <iostream>

using namespace std;
using namespace Victor;
using namespace Victor::Biopool;

/*
 * 
 */
int main(int argc, char** argv) {
    //Rotator* rotator = new KabschMethod;
    //Load proteins
    //string proteine1 = "/home/cecco/Desktop/TestBIO2/15C8_H_input.pdb";
    //string proteine2 = "/home/cecco/Desktop/TestBIO2/25C8_H_input.pdb";
    string proteine1 = "/home/cecco/Desktop/AVANZATITestBIO2/T0760TS008_1";
    string proteine2 = "/home/cecco/Desktop/AVANZATITestBIO2/T0760TS008_2";
    ifstream inFile1(proteine1.c_str());
    ifstream inFile2(proteine2.c_str());
    // creates the PdbLoader1 object
    PdbLoader pl1(inFile1);
    // creates the PdbLoader2 object
    PdbLoader pl2(inFile2);


    Protein* prot1 = new Protein();
    Protein* prot2 = new Protein();

    prot1->load(pl1);
    prot2->load(pl2);

    SuperImpositor* superImpositor = new SuperImpositor(prot1, prot2, "");
    double rmsd = superImpositor->calculateRMSD();
    double maxSub = superImpositor->calculateMaxSub();
    cout << "l'rmsd e':" << rmsd << "\n";
    cout << "il maxsub e':" << maxSub << "\n";
    Spacer newSet1 = superImpositor->getRMSDset1();
    Spacer newSet2 = superImpositor->getRMSDset2();
    string proteineOUTPUT1 = "/home/cecco/Desktop/TestBIO2/output1.pdb";
    string proteineOUTPUT2 = "/home/cecco/Desktop/TestBIO2/output2.pdb";
    ofstream outFile1(proteineOUTPUT1.c_str());
    ofstream outFile2(proteineOUTPUT2.c_str());

    PdbSaver saveSet1(outFile1);
    PdbSaver saveSet2(outFile2);
    saveSet1.saveSpacer(newSet1);
    saveSet2.saveSpacer(newSet2);
    return 0;
}

