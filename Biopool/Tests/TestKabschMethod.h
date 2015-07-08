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

#include <KabschMethod.h>

#include <PdbLoader.h>

using namespace std;
using namespace Victor::Biopool;

class TestKabschMethod : public CppUnit::TestFixture {
private:
    Spacer *testKabschMethod;
public:

    TestKabschMethod() : testKabschMethod(NULL) {
    }

    virtual ~TestKabschMethod() {
        delete testKabschMethod;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestKabschMethod");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestKabschMethod>("Test1 - calculate rototraslation matrix.",
                &TestKabschMethod::testTestKabschMethod_A));

        return suiteOfTests;
    }

    /// Setup method

    void setUp() {




    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testTestKabschMethod_A() {
        Rotator* rotationAlgorith = new KabschMethod();
        Eigen::Matrix3Xd set1(3, 2), set2(3, 2);

        set1(0, 0) = 1;
        set1(0, 1) = 2;
        set1(1, 0) = 0;
        set1(1, 1) = 0;
        set1(2, 0) = 0;
        set1(2, 1) = 0;

        set2(0, 0) = 0;
        set2(0, 1) = -1;
        set2(1, 0) = 1;
        set2(1, 1) = 1;
        set2(2, 0) = 1;
        set2(2, 1) = 1;

        Eigen::Affine3d* input = rotationAlgorith->rotate(set1, set2);
        cout << "\n\nThe rotation matrix is:\n";
        Eigen::Matrix3d rot = input->linear();
        cout << rot << "\n";
        cout << "The translation vector is:\n";
        Eigen::Vector3d tra = input->translation();
        cout << tra << "\n";
        CPPUNIT_ASSERT(rot(0, 0) == -1 && rot(1, 1) == 1 && rot(2, 2) == -1 && tra(0) == 1 && tra(1) == 1 && tra(2) == 1);
    }

};