/* 
 * File:   superImpositor.cc
 * Author: cecco
 * 
 * Created on June 8, 2015, 4:22 PM
 */

#include "SuperImpositor.h"
#include <math.h>
#include <vector>

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
    matrixSet1 = fromSpacerToMatrix3Xd(*(set1));
    matrixSet2 = fromSpacerToMatrix3Xd(*(set2));
}

SuperImpositor::~SuperImpositor() {
}

double SuperImpositor::calculateRMSD() {
    Eigen::Affine3d* rotoTraslation = rotationAlgorith->rotate(matrixSet1, matrixSet2);
    calculateRotation(matrixSet1, matrixSet2, rotoTraslation);
    RMSDset1 = fromMatrix3XdToSpacer(matrixSet1, 1);
    RMSDset2 = fromMatrix3XdToSpacer(matrixSet2, 2);
    //Calculate rmsd
    double distanceSquared = 0;
    double msd;
    int columns = matrixSet1.cols();
    for (int col = 0; col < columns; col++) {
        //distanceSquared = distanceSquared + pow((matrixSet1(0, col) - matrixSet2(0, col)), 2);
        //distanceSquared = distanceSquared + pow((matrixSet1(1, col) - matrixSet2(1, col)), 2);
        //distanceSquared = distanceSquared + pow((matrixSet1(2, col) - matrixSet2(2, col)), 2);
        //distance=sqrt(distanceSquared);
        distanceSquared += (matrixSet1.col(col) - matrixSet2.col(col)).squaredNorm();

        //rmsdSQ+=distance^2;
        //distance=0;
    }
    msd = distanceSquared / columns;

    return sqrt(msd);
}

double SuperImpositor::calculateMaxSub() {
    std::vector<std::pair<int, int> > vectorSet;
    for (int i = 0; i < matrixSet1.cols(); i++) {
        vectorSet.push_back(std::make_pair(i + 1, i + 1));
    }
    //M is a general vector. Is possible to use a general alignment  and not only (x,x),(y,y),ecc...
    maxSubAlignment(matrixSet1, matrixSet2, vectorSet);

    return 0;
}

std::vector<std::pair<int, int> > SuperImpositor::maxSubAlignment(Eigen::Matrix3Xd& firstSet, Eigen::Matrix3Xd& secondSet, std::vector<std::pair<int, int> > vectorSet) {
    double distance = 3.5;
    long unsigned int sMax = 0;
    int n = firstSet.cols();
    int L = 4;
    std::vector<std::pair<int, int> > mMax;
    std::vector<std::pair<int, int> > M;
    //nb ricorda che non e' detto che siano consecutivi gli elementi mentre estraggo nei vari cicli di Extend
    //std::vector<std::pair<int, int> >::iterator it;
    //it = vectorSet.begin();
    for (int i = 0; i < n - L + 1; i++) {
        //Add the L elements
        for (int i = 0; i < L; i++) {
            M.push_back(vectorSet[i]);
        }
        M = Extend(M, firstSet, secondSet, distance, L);
        if (M.size() > sMax) {

            sMax = M.size();
            mMax = M;
        }
    }
    return mMax;
}

std::vector< std::pair<int, int> > SuperImpositor::Extend(std::vector<std::pair<int, int> > M, Eigen::Matrix3Xd A, Eigen::Matrix3Xd B, double d, int L) {
    std::vector<std::pair<int, int> > n;
    int k = 4;
    //    Eigen::Matrix3Xd M1(3, L);
    //    Eigen::Matrix3Xd M2(3, L);
    Eigen::Matrix3Xd N1;
    Eigen::Matrix3Xd N2;
    for (int i = 0; i < L; i++) {
        N1 << N1, A.col(M[i].first);
        N2 << N2, B.col(M[i].second);
    }
    
    
    //BISOGNA TENERE TRACCIA DEGLI ELEMENTI PRESENTI IN N1 perche' dopo un ciclo A puo' non avvere un elemento .col(M[i].first)
    
    
    
    for (int j = 1; j <= k; j++) {
        Eigen::Affine3d* rotoTraslation = rotationAlgorith->rotate(N1, N2);
        calculateRotation(A, B, rotoTraslation);

        double distance;
        //Calculate the distance between the points

        for (int i = 0; i < A.cols(); i++) {
            distance = (A.col(M[i].first) - B.col(M[i].second)).squaredNorm();
            //            distance = pow((A(0, i) - B(0, i)), 2);
            //            distance = distance + pow((A(1, i) - B(1, i)), 2);
            //            distance = distance + pow((A(2, i) - B(2, i)), 2);
            distance = sqrt(distance);
            if (distance < ((j * d) / k)) {
                N1 << A.col(M[i].first);
                N2 << B.col(M[i].second);
            }
        }
    }
    Eigen::Affine3d* rotoTraslation = rotationAlgorith->rotate(N1, N2);
    calculateRotation(A, B, rotoTraslation);
    double distance;
    //Calculate the distance between the points

    for (int i = 0; i < A.cols(); i++) {
        distance = (A.col(M[i].first) - B.col(M[i].second)).squaredNorm();
        //            distance = pow((A(0, i) - B(0, i)), 2);
        //            distance = distance + pow((A(1, i) - B(1, i)), 2);
        //            distance = distance + pow((A(2, i) - B(2, i)), 2);
        distance = sqrt(distance);
        if (distance < d) {
            N1 << A.col(M[i].first);
            N2 << B.col(M[i].second);
        }
    }

    return M;
}

void SuperImpositor::calculateRotation(Eigen::Matrix3Xd& firstSet, Eigen::Matrix3Xd& secondSet, Eigen::Affine3d* rotoTraslation) {
    Eigen::Matrix3Xd R = rotoTraslation->linear();
    Eigen::Vector3d S = rotoTraslation->translation();
    //Apply the rototraslation
    int NumAmino = firstSet.cols();
    Eigen::Matrix3Xd rotoTraslSet(3, NumAmino);
    for (int col = 0; col < NumAmino; col++) {

        rotoTraslSet.col(col) = R * firstSet.col(col) + S;
    }
    //Save rototraslation
    firstSet = rotoTraslSet;
}

void SuperImpositor::rotationApplication(Eigen::Matrix3Xd R) {

    vgMatrix3<double> rotationMatrix;
    rotationMatrix = fromMatrix3XdTovgMatrix3(R);

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
    }

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

Spacer SuperImpositor::fromMatrix3XdToSpacer(Eigen::Matrix3Xd matrix3Xd, int num) const {
    //Change atoms coordinates of set1
    Spacer newSpacer;
    if (num == 1)
        newSpacer = *(set1);
    else
        newSpacer = *(set2);
    int NumAmino = matrix3Xd.cols();
    for (int i = 0; i < NumAmino; i++) {
        newSpacer.getAmino(i)[CA].setCoords(fromVector3dTovgVector3(matrix3Xd.col(i)));
    }
    return newSpacer;
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

Spacer SuperImpositor::getRMSDset2() const {
    return RMSDset2;
}

Spacer SuperImpositor::getRMSDset1() const {
    return RMSDset1;
}