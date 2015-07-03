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
#include <GetArg.h>

using namespace std;
using namespace Victor;
using namespace Victor::Biopool;



//TODO:
//-CONTROLLARE NELLE REGOLE DEL PROGETTO COME LEGGERE GLI INPUT
//-CONTROLLARE SE VUOLE SEMPRE TUTTI GLI OUTPUT O SOLO UN METODO RICHIESTO
//-CHIEDERE SE VA BENE L' APPROSSIMAZIONE OTTENUTA DI 10^-2 PER GDT E TMSCORE o chiedere come si implementa (tipo per l'erase)+ chiedi se basta test per le 4 funz principali
//-STAMPARE IL PDF PER OGNI METODO, prima + seconda proteina una modificata
//-SCRIVERE la documentazione
//-SCRIVERE TEST
//

void sShowHelp() {
    cout << "\n"
            << "   Options: \n"
            << "\t-i <filename>      \t\t Input file for PDB structure \n"
            << "\n"

            << "\t-r rmsd            \t\t run the rmsd algorithm and display the value \n"
            << "\t-m maxsub          \t\t run the maxsub algorithm and display the value \n"
            << "\t-g gdt             \t\t run the gdt algorithm and display the value \n"
            << "\t-t tmscore         \t\t run the tmscore algorithm and display the value \n"
            << "\n"

            << "\t-R rmsd            \t\t get in the output file the rmsd rototraslated proteins \n"
            << "\t-M maxsub          \t\t get in the output file the maxsub rototraslated proteins \n"
            << "\t-G gdt             \t\t get in the output file the gdt rototraslated proteins \n"
            << "\t-T tmscore         \t\t get in the output file the tmscore rototraslated proteins \n"
            << "\n"

            << "\t-x rotation method \t\t select the rotation method from: \n"
            << "\t                   \t\t -kabsch method \n"

            << "\n";
}

void savePdbOutput(vector <Spacer> spacers, string name);

int main(int nArgs, char* argv[]) {

    if (getArg("h", nArgs, argv)) {
        sShowHelp();
        return 1;
    };

    vector<string> inputFile;
    getArg("i", inputFile, nArgs, argv, "!");
    if (inputFile.size() != 2) {
        cout << "Error you need to get in input exactly 2 input file. Aborting. (-h for help)" << endl;
        return -1;
    }

    bool rmsd, maxsub, gdt, tmscore;
    bool rmsdOutput, maxsubOutput, gdtOutput, tmscoreOutput;

    rmsd = getArg("r", nArgs, argv);
    maxsub = getArg("m", nArgs, argv);
    gdt = getArg("g", nArgs, argv);
    tmscore = getArg("t", nArgs, argv);

    rmsdOutput = getArg("R", nArgs, argv);
    maxsubOutput = getArg("M", nArgs, argv);
    gdtOutput = getArg("G", nArgs, argv);
    tmscoreOutput = getArg("T", nArgs, argv);

    if (rmsd || maxsub || gdt || tmscore || rmsdOutput || maxsubOutput || gdtOutput || tmscoreOutput) {
        string rotationMethod;
        getArg("x", rotationMethod, nArgs, argv, "");

        ifstream inFile1(inputFile[0].c_str());
        ifstream inFile2(inputFile[1].c_str());

        // creates the PdbLoader objects
        PdbLoader pl1(inFile1);
        PdbLoader pl2(inFile2);

        Protein* prot1 = new Protein();
        Protein* prot2 = new Protein();

        prot1->load(pl1);
        prot2->load(pl2);

        SuperImpositor* superImpositor = new SuperImpositor(prot1, prot2, rotationMethod);

        if (rmsd) {
            double rmsd = superImpositor->calculateRMSD();

            if (!(rmsdOutput || maxsubOutput || gdtOutput || tmscoreOutput)) {
                vector <Spacer> spacers;
                Spacer newSet1 = superImpositor->getRMSDset1();
                Spacer newSet2 = superImpositor->getRMSDset2();

                spacers.push_back(newSet1);
                spacers.push_back(newSet2);
                savePdbOutput(spacers, "rmsd");
            }
            cout << "The rmsd value is': " << rmsd << "\n";
        }

        if (rmsdOutput) {
            superImpositor->calculateRMSD();
            vector <Spacer> spacers;
            Spacer newSet1 = superImpositor->getRMSDset1();
            Spacer newSet2 = superImpositor->getRMSDset2();

            spacers.push_back(newSet1);
            spacers.push_back(newSet2);
            savePdbOutput(spacers, "rmsd");
        }

        if (maxsub) {
            std::vector<std::pair<int, int> > vectorSet;

            for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
                vectorSet.push_back(std::make_pair(i, i));
            }
            double maxSub = superImpositor->calculateMaxSub(3.5, vectorSet, 'y');

            if (!(rmsdOutput || maxsubOutput || gdtOutput || tmscoreOutput)) {
                vector <Spacer> spacers;
                Spacer newSet1 = superImpositor->getMaxSubset1();
                Spacer newSet2 = superImpositor->getMaxSubset2();

                spacers.push_back(newSet1);
                spacers.push_back(newSet2);
                savePdbOutput(spacers, "maxsub");
            }
            cout << "The maxsub value is: " << maxSub << "\n";
        }

        if (maxsubOutput) {
            std::vector<std::pair<int, int> > vectorSet;

            for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
                vectorSet.push_back(std::make_pair(i, i));
            }
            superImpositor->calculateMaxSub(3.5, vectorSet, 'y');

            vector <Spacer> spacers;
            Spacer newSet1 = superImpositor->getMaxSubset1();
            Spacer newSet2 = superImpositor->getMaxSubset2();
            spacers.push_back(newSet1);
            spacers.push_back(newSet2);
            savePdbOutput(spacers, "maxsub");
        }

        if (gdt) {
            std::vector<std::pair<int, int> > vectorSet;

            for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
                vectorSet.push_back(std::make_pair(i, i));
            }
            double gdt = superImpositor->calculateGdt(vectorSet);
            if (!(rmsdOutput || maxsubOutput || gdtOutput || tmscoreOutput)) {
                vector <Spacer> spacers;
                Spacer newSet1 = superImpositor->getGdtset1_1();
                Spacer newSet2 = superImpositor->getGdtset1_2();

                Spacer newSet3 = superImpositor->getGdtset2_1();
                Spacer newSet4 = superImpositor->getGdtset2_2();

                Spacer newSet5 = superImpositor->getGdtset3_1();
                Spacer newSet6 = superImpositor->getGdtset3_2();

                Spacer newSet7 = superImpositor->getGdtset4_1();
                Spacer newSet8 = superImpositor->getGdtset4_2();

                spacers.push_back(newSet1);
                spacers.push_back(newSet2);
                spacers.push_back(newSet3);
                spacers.push_back(newSet4);
                spacers.push_back(newSet5);
                spacers.push_back(newSet6);
                spacers.push_back(newSet7);
                spacers.push_back(newSet8);
                //savePdbOutput(spacers, "gdt");
            }
            cout << "The gdt value is: " << gdt << "\n";
        }

        if (gdtOutput) {
            std::vector<std::pair<int, int> > vectorSet;

            for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
                vectorSet.push_back(std::make_pair(i, i));
            }
            superImpositor->calculateGdt(vectorSet);

            vector <Spacer> spacers;
            Spacer newSet1 = superImpositor->getGdtset1_1();
            Spacer newSet2 = superImpositor->getGdtset1_2();

            Spacer newSet3 = superImpositor->getGdtset2_1();
            Spacer newSet4 = superImpositor->getGdtset2_2();

            Spacer newSet5 = superImpositor->getGdtset3_1();
            Spacer newSet6 = superImpositor->getGdtset3_2();

            Spacer newSet7 = superImpositor->getGdtset4_1();
            Spacer newSet8 = superImpositor->getGdtset4_2();

            spacers.push_back(newSet1);
            spacers.push_back(newSet2);
            spacers.push_back(newSet3);
            spacers.push_back(newSet4);
            spacers.push_back(newSet5);
            spacers.push_back(newSet6);
            spacers.push_back(newSet7);
            spacers.push_back(newSet8);
            savePdbOutput(spacers, "gdt");
        }

        if (tmscore) {
            std::vector<std::pair<int, int> > vectorSet;

            for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
                vectorSet.push_back(std::make_pair(i, i));
            }
            double TMScore = superImpositor->calculateTMScore(vectorSet);
            if (!(rmsdOutput || maxsubOutput || gdtOutput || tmscoreOutput)) {
                vector <Spacer> spacers;
                Spacer newSet1 = superImpositor->getTMScoreset1();
                Spacer newSet2 = superImpositor->getTMScoreset2();

                spacers.push_back(newSet1);
                spacers.push_back(newSet2);
                savePdbOutput(spacers, "TMScore");
            }
            cout << "The TMScore value is: " << TMScore << "\n";
        }

        if (tmscoreOutput) {
            std::vector<std::pair<int, int> > vectorSet;

            for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
                vectorSet.push_back(std::make_pair(i, i));
            }
            superImpositor->calculateTMScore(vectorSet);

            vector <Spacer> spacers;
            Spacer newSet1 = superImpositor->getTMScoreset1();
            Spacer newSet2 = superImpositor->getTMScoreset2();

            spacers.push_back(newSet1);
            spacers.push_back(newSet2);
            savePdbOutput(spacers, "TMScore");
        }


    } else {
        cout << "No method select (-h for help)\n";
    }

    return 0;
}

void savePdbOutput(vector <Spacer> spacers, string name) {
    string proteineOUTPUT = name + "output.pdb";
    ofstream outFile(proteineOUTPUT.c_str());

    PdbSaver saveSet(outFile);
    for (unsigned int i = 0; i < spacers.size(); i++) {
        saveSet.saveSpacer(spacers[i]);
    }
}
