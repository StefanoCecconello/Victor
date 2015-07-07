/*
 * TestAtom.cpp
 *
 *  Created on: Oct 6th, 2014
 *      Author: Layla Hirsh
 */

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <SuperImpositor.h>

#include <PdbLoader.h>

using namespace std;
using namespace Victor::Biopool;

class TestSuperImpositor : public CppUnit::TestFixture {
private:
    Spacer *testSuperImpositor;
public:

    TestSuperImpositor() : testSuperImpositor(NULL) {
    }

    virtual ~TestSuperImpositor() {
        delete testSuperImpositor;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestSuperImpositor");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestSuperImpositor>("Test1 - calculate rmsd.",
                &TestSuperImpositor::testTestSuperImpositor_A));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestSuperImpositor>("Test2 - calculate maxsub.",
                &TestSuperImpositor::testTestSuperImpositor_B));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestSuperImpositor>("Test2 - calculate maxsub.",
                &TestSuperImpositor::testTestSuperImpositor_C));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestSuperImpositor>("Test2 - calculate maxsub.",
                &TestSuperImpositor::testTestSuperImpositor_D));


        return suiteOfTests;
    }

    /// Setup method

    void setUp() {




    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testTestSuperImpositor_A() {
        string path = getenv("VICTOR_ROOT");
        string inputFile1 = path + "Biopool/Tests/data/TestAlign1.pdb";
        string inputFile2 = path + "Biopool/Tests/data/TestAlign2.pdb";

        ifstream inFile1(inputFile1.c_str());
        ifstream inFile2(inputFile2.c_str());
        if (!inFile1)
            ERROR("First file not found.", exception);
        if (!inFile2)
            ERROR("Second File not found.", exception);
        PdbLoader pl1(inFile1);
        PdbLoader pl2(inFile2);
        Protein* prot1 = new Protein();
        Protein* prot2 = new Protein();
        pl1.setNoVerbose();
        pl2.setNoVerbose();
        pl1.setNoHAtoms();
        pl2.setNoHAtoms();
        prot1->load(pl1);
        prot2->load(pl2);
        string method = "Kabsch";
        SuperImpositor* superImpositor = new SuperImpositor(prot1, prot2, method);
        superImpositor->calculateRMSD();
        double rmsd = superImpositor->getRmsdValue();
        cout << "\n\nThe rmsd value between the two protein is: " << rmsd << "\n";
        CPPUNIT_ASSERT((rmsd - 21.732) < EPSILON && (-(rmsd - 21.732)) < EPSILON);
    }

    void testTestSuperImpositor_B() {
        string path = getenv("VICTOR_ROOT");
        string inputFile1 = path + "Biopool/Tests/data/TestAlign1.pdb";
        string inputFile2 = path + "Biopool/Tests/data/TestAlign2.pdb";

        ifstream inFile1(inputFile1.c_str());
        ifstream inFile2(inputFile2.c_str());
        if (!inFile1)
            ERROR("First file not found.", exception);
        if (!inFile2)
            ERROR("Second File not found.", exception);
        PdbLoader pl1(inFile1);
        PdbLoader pl2(inFile2);
        Protein* prot1 = new Protein();
        Protein* prot2 = new Protein();
        pl1.setNoVerbose();
        pl2.setNoVerbose();
        pl1.setNoHAtoms();
        pl2.setNoHAtoms();
        prot1->load(pl1);
        prot2->load(pl2);
        string method = "Kabsch";

        SuperImpositor* superImpositor = new SuperImpositor(prot1, prot2, method);

        std::vector<std::pair<int, int> > vectorSet;

        for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
            vectorSet.push_back(std::make_pair(i, i));
        }

        superImpositor->calculateMaxSub(3.5, vectorSet, 'y');
        double maxSub = superImpositor->getMaxsubValue();
        cout << "\n\nThe maxsub value between the two protein is: " << maxSub << "\n";
        CPPUNIT_ASSERT((maxSub - 0.626291) < EPSILON && (-(maxSub - 0.626291)) < EPSILON);
    }

    void testTestSuperImpositor_C() {
        string path = getenv("VICTOR_ROOT");
        string inputFile1 = path + "Biopool/Tests/data/TestAlign1.pdb";
        string inputFile2 = path + "Biopool/Tests/data/TestAlign2.pdb";

        ifstream inFile1(inputFile1.c_str());
        ifstream inFile2(inputFile2.c_str());
        if (!inFile1)
            ERROR("First file not found.", exception);
        if (!inFile2)
            ERROR("Second File not found.", exception);
        PdbLoader pl1(inFile1);
        PdbLoader pl2(inFile2);
        Protein* prot1 = new Protein();
        Protein* prot2 = new Protein();
        pl1.setNoVerbose();
        pl2.setNoVerbose();
        pl1.setNoHAtoms();
        pl2.setNoHAtoms();
        prot1->load(pl1);
        prot2->load(pl2);
        string method = "Kabsch";

        SuperImpositor* superImpositor = new SuperImpositor(prot1, prot2, method);

        std::vector<std::pair<int, int> > vectorSet;

        for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
            vectorSet.push_back(std::make_pair(i, i));
        }

        superImpositor->calculateGdt(vectorSet);
        double gdt = superImpositor->getGdtValue();
        cout << "\n\nThe gdt value between the two protein is: " << gdt << "\n";
        CPPUNIT_ASSERT(((gdt < 0.6725) + 0.1) && (gdt > (0.6725 - 0.1)));
    }

    void testTestSuperImpositor_D() {
        string path = getenv("VICTOR_ROOT");
        string inputFile1 = path + "Biopool/Tests/data/TestAlign1.pdb";
        string inputFile2 = path + "Biopool/Tests/data/TestAlign2.pdb";

        ifstream inFile1(inputFile1.c_str());
        ifstream inFile2(inputFile2.c_str());
        if (!inFile1)
            ERROR("First file not found.", exception);
        if (!inFile2)
            ERROR("Second File not found.", exception);
        PdbLoader pl1(inFile1);
        PdbLoader pl2(inFile2);
        Protein* prot1 = new Protein();
        Protein* prot2 = new Protein();
        pl1.setNoVerbose();
        pl2.setNoVerbose();
        pl1.setNoHAtoms();
        pl2.setNoHAtoms();
        prot1->load(pl1);
        prot2->load(pl2);
        string method = "Kabsch";

        SuperImpositor* superImpositor = new SuperImpositor(prot1, prot2, method);

        std::vector<std::pair<int, int> > vectorSet;

        for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
            vectorSet.push_back(std::make_pair(i, i));
        }

        superImpositor->calculateTMScore(vectorSet);
        double TMScore = superImpositor->getTMScoreValue();
        cout << "\n\nThe TMScore value between the two protein is: " << TMScore << "\n";
        CPPUNIT_ASSERT(((TMScore < 0.7408) + 0.1) && (TMScore > (0.7408 - 0.1)));
    }


};