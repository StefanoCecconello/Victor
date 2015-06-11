/* 
 * File:   main.cc
 * Author: cecco
 *
 * Created on May 29, 2015, 4:18 PM
 */

#include <cstdlib>
#include <SuperImpositor.h>
#include <PdbLoader.h>
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
    string proteine1 = "/home/cecco/Desktop/TestBIO2/15C8_H_input.pdb";
    string proteine2 = "/home/cecco/Desktop/TestBIO2/25C8_H_input.pdb";
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

    SuperImpositor* aaa = new SuperImpositor(prot1, prot2, "");

    return 51;
}

