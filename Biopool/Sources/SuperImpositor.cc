/* 
 * File:   superImpositor.cc
 * Author: cecco
 * 
 * Created on June 8, 2015, 4:22 PM
 */

#include "SuperImpositor.h"


using namespace Victor;
using namespace Victor::Biopool;

SuperImpositor::SuperImpositor(Protein* firstProtein, Protein* secondProtein, string method = "") {
    if (method == "Kabsch") {
        rotationAlgorith = new KabschMethod();
    }

    //set default rotation method
    if (method == "") {
        rotationAlgorith = new KabschMethod();
    }

    //Get the spacers
    set1 = firstProtein->getSpacer((unsigned int) 0);
    set2 = secondProtein->getSpacer((unsigned int) 0);
}

SuperImpositor::~SuperImpositor() {
}

double SuperImpositor::calculateRMSD() {
    int NumAmino = set1->sizeAmino();


    Eigen::Matrix3Xd matrixSet1 = fromSpacerToMatrix3Xd(*(set1));
    Eigen::Matrix3Xd matrixSet2 = fromSpacerToMatrix3Xd(*(set2));
    Eigen::Affine3d* rotoTraslation = rotationAlgorith->rotate(matrixSet1, matrixSet2);

    Eigen::Matrix3Xd R = rotoTraslation->linear();
    Eigen::Vector3d S = rotoTraslation->translation();


    //translationApplication(S);

    //Eigen::Matrix3Xd a;
    //a=R*fromSpacerToMatrix3Xd(set1);
    Eigen::Matrix3Xd rotoTraslSet(3, NumAmino);
    for (int col = 0; col < matrixSet1.cols(); col++)
        rotoTraslSet.col(col) = R * matrixSet1.col(col) + S;

    //    cout << R << "\n" << "\n" << "\n" << "\n";
    //    cout << S << "\n" << "\n" << "\n" << "\n";
    //
    //
    //    cout << fromSpacerToMatrix3Xd(set1) << "\n" << "\n" << "\n" << "\n";
    //    cout << fromSpacerToMatrix3Xd(set2) << "\n" << "\n" << "\n" << "\n";

    for (int i = 0; i < NumAmino; i++) {
        set2->getAmino(i)[CA].setCoords(fromVector3dTovgVector3(rotoTraslSet.col(i)));
    }
    fromSpacerToMatrix3Xd(*(set1));
    fromSpacerToMatrix3Xd(*(set2));

    //for (int col = 0; col < matrixSet1.cols(); col++)
    //    matrixSet2.col(col) = R * matrixSet1.col(col) + S;


    //cout << fromSpacerToMatrix3Xd(set2) << "\n" << "\n" << "\n" << "\n";
    //set2->addRot(rotationMatrix);
    //set2->addTrans(traslationVector);

    //rotationApplication(R);
    //translationApplication(S);
    //cout << fromSpacerToMatrix3Xd(set1) << "\n" << "\n" << "\n" << "\n";
    //cout << fromSpacerToMatrix3Xd(set2) << "\n" << "\n" << "\n" << "\n";

    return 0;
}

void SuperImpositor::rotationApplication(Eigen::Matrix3Xd R) {
    vgMatrix3<double> rotationMatrix;
    rotationMatrix = fromMatrix3XdTovgMatrix3(R);


    //    int NumAmino = set2->sizeAmino();
    //    for (int i = 0; i < NumAmino; i++) {
    //        set2->getAmino(i).addRot(rotationMatrix);
    //    }

}

void SuperImpositor::translationApplication(Eigen::Vector3d S) {
    vgVector3<double> traslationVector;
    traslationVector = fromVector3dTovgVector3(S);

    int NumAmino = set2->sizeAmino();
    for (int i = 0; i < NumAmino; i++) {
        set2->getAmino(i).addTrans(traslationVector);
    }
}

Eigen::Matrix3Xd SuperImpositor::fromSpacerToMatrix3Xd(Spacer spacerSet) const {
    int NumAmino = spacerSet.sizeAmino();
    Atom CAAtoms[NumAmino];
    Eigen::Matrix3Xd matrixSet(3, NumAmino);

    for (int i = 0; i < NumAmino; i++) {
        CAAtoms[i] = spacerSet.getAmino(i)[CA];
    }

    for (int col = 0; col < NumAmino; col++) {
        vgVector3<double> coords = CAAtoms[col].getCoords();
        matrixSet(0, col) = coords[0];
        matrixSet(1, col) = coords[1];
        matrixSet(2, col) = coords[2];
        cout << coords[0] << "\t";
        cout << coords[1] << "\t";
        cout << coords[2] << "\n";
    }
    cout << "\n" << "\n" << "\n" << "\n";

    return matrixSet;
}

vgMatrix3<double> SuperImpositor::fromMatrix3XdTovgMatrix3(Eigen::Matrix3Xd matrix3Xd) const {
    //Conversion from Matrix3Xd to vgMatrix3
    vgMatrix3<double> rotationMatrix;
    rotationMatrix.x.x = matrix3Xd(0, 0);
    rotationMatrix.x.y = matrix3Xd(0, 1);
    rotationMatrix.x.z = matrix3Xd(0, 2);
    rotationMatrix.y.x = matrix3Xd(1, 0);
    rotationMatrix.y.y = matrix3Xd(1, 1);
    rotationMatrix.y.z = matrix3Xd(1, 2);
    rotationMatrix.z.x = matrix3Xd(2, 0);
    rotationMatrix.z.y = matrix3Xd(2, 1);
    rotationMatrix.z.z = matrix3Xd(2, 2);
    return rotationMatrix;
}

vgVector3<double> SuperImpositor::fromVector3dTovgVector3(Eigen::Vector3d Vector3d) const {
    vgVector3<double> traslationVector;
    for (int r = 0; r < 3; r++) {
        traslationVector[r] = Vector3d(r);
    }
    return traslationVector;
}

Spacer* SuperImpositor::getSet2() const {
    return set2;
}

Spacer* SuperImpositor::getSet1() const {
    return set1;
}