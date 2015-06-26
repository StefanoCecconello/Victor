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
    Eigen::Matrix3Xd modifyMatrixSet1 = matrixSet1;
    Eigen::Matrix3Xd modifyMatrixSet2 = matrixSet2;
    Eigen::Affine3d* rotoTraslation = rotationAlgorith->rotate(modifyMatrixSet1, modifyMatrixSet2);
    calculateRotation(modifyMatrixSet1, modifyMatrixSet2, rotoTraslation);
    RMSDset1 = fromMatrix3XdToSpacer(modifyMatrixSet1, 1);
    RMSDset2 = fromMatrix3XdToSpacer(modifyMatrixSet2, 2);
    //Calculate rmsd
    double distanceSquared = 0;
    double msd;
    int columns = modifyMatrixSet1.cols();
    for (int col = 0; col < columns; col++) {
        //distanceSquared = distanceSquared + pow((modifyMatrixSet1(0, col) - modifyMatrixSet2(0, col)), 2);
        //distanceSquared = distanceSquared + pow((modifyMatrixSet1(1, col) - modifyMatrixSet2(1, col)), 2);
        //distanceSquared = distanceSquared + pow((modifyMatrixSet1(2, col) - modifyMatrixSet2(2, col)), 2);
        //distanceSquared=sqrt(distanceSquared);
        distanceSquared += (modifyMatrixSet1.col(col) - modifyMatrixSet2.col(col)).squaredNorm();

        //rmsdSQ+=distanceSquared^2;
        //distance=0;
    }
    msd = distanceSquared / columns;

    return sqrt(msd);
}

double SuperImpositor::calculateMaxSub(double d, std::vector<std::pair<int, int> > vectorSet, char E) {



    double distance;
    double sum = 0;

    //M is a general vector. Is possible to use a general alignment  and not only (x,x),(y,y),ecc...
    std::vector<std::pair<int, int> > mMax;
    Eigen::Matrix3Xd modifyMatrixSet1 = matrixSet1;
    Eigen::Matrix3Xd modifyMatrixSet2 = matrixSet2;
    mMax = maxSubAlignment(modifyMatrixSet1, modifyMatrixSet2, vectorSet, d, E);
    for (unsigned int i = 0; i < mMax.size(); i++) {
        //distance = (matrixSet1.col(mMax[i].first) - matrixSet2.col(mMax[i].second)).squaredNorm();
        distance = (modifyMatrixSet1.col(mMax[i].first) - modifyMatrixSet2.col(mMax[i].second)).squaredNorm();
        distance = sqrt(distance);
        sum = sum + (1 / (1 + pow((distance / d), 2)));
    }
    return sum / modifyMatrixSet1.cols();
}

double SuperImpositor::calculateGdt(std::vector<std::pair<int, int> > vectorSet) {
    //0.6135 con maxsub programmino per 4 e 5
    //    cout<<calculateMaxSub(1, vectorSet)<<"\n";
    //    cout<<calculateMaxSub(2, vectorSet)<<"\n";
    //    cout<<calculateMaxSub(4, vectorSet)<<"\n";
    //    cout<<calculateMaxSub(8, vectorSet)<<"\n";
    return (calculateMaxSub(1, vectorSet, 'n') + calculateMaxSub(2, vectorSet, 'n') + calculateMaxSub(4, vectorSet, 'n') + calculateMaxSub(8, vectorSet, 'n')) / 4;
}

double SuperImpositor::calculateTMScore(std::vector<std::pair<int, int> > vectorSet) {
    Eigen::Matrix3Xd modifyMatrixSet1 = matrixSet1;
    //    Eigen::Matrix3Xd modifyMatrixSet2 = matrixSet2;
    //    Eigen::Affine3d* rotoTraslation = rotationAlgorith->rotate(modifyMatrixSet1, modifyMatrixSet2);
    //    calculateRotation(modifyMatrixSet1, modifyMatrixSet2, rotoTraslation);
    //    double sum = 0;
    //    double distance;
    double d0;
    d0 = (1.24 * cbrt(modifyMatrixSet1.cols() - 15)) - 1.8;
    //    for (unsigned int i = 0; i < modifyMatrixSet2.cols(); i++) {
    //        distance = (modifyMatrixSet1.col(vectorSet[i].first) - modifyMatrixSet2.col(vectorSet[i].second)).squaredNorm();
    //        distance = sqrt(distance);
    //        sum = sum + (1 / (1 + pow((distance / d0), 2)));
    //    }
    //    sum = sum / modifyMatrixSet1.cols();

    return calculateMaxSub(d0, vectorSet, 'n');
}

std::vector<std::pair<int, int> > SuperImpositor::maxSubAlignment(Eigen::Matrix3Xd& firstSet, Eigen::Matrix3Xd& secondSet, std::vector<std::pair<int, int> > vectorSet, double d, char E) {
    long unsigned int sMax = 0;
    int n = firstSet.cols();
    int L = 4;
    Eigen::Matrix3Xd modifiedFirstSet;
    Eigen::Matrix3Xd modifiedSecondSet;
    Eigen::Matrix3Xd optFirstSet;
    Eigen::Matrix3Xd optSecondSet;
    std::vector<std::pair<int, int> > mMax;
    std::vector<std::pair<int, int> > M;
    for (int i = 0; i < n - L + 1; i++) {
        modifiedFirstSet = firstSet;
        modifiedSecondSet = secondSet;
        M.clear();
        //Add the L elements
        for (int j = 0; j < L; j++) {
            M.push_back(vectorSet[i + j]);
        }
        M = Extend(M, vectorSet, modifiedFirstSet, modifiedSecondSet, d, L, n, E);
        if (M.size() > sMax) {
            optFirstSet = modifiedFirstSet;
            optSecondSet = modifiedSecondSet;
            sMax = M.size();
            mMax = M;
        }
    }
    firstSet = optFirstSet;
    secondSet = optSecondSet;
    return mMax;
}

std::vector< std::pair<int, int> > SuperImpositor::Extend(std::vector<std::pair<int, int> > M, std::vector<std::pair<int, int> > vectorSet, Eigen::Matrix3Xd& A, Eigen::Matrix3Xd& B, double d, int L, int n, char E) {
    std::vector<std::pair<int, int> > N;
    double threshold;
    int k = 4;
    double distance;
    Eigen::Matrix3Xd M1(3, L);
    Eigen::Matrix3Xd M2(3, L);

    for (int i = 0; i < L; i++) {
        M1.col(i) = A.col(M[i].first);
        M2.col(i) = B.col(M[i].second);
    }

    Eigen::Matrix3Xd ARototrasled;
    Eigen::Matrix3Xd BRototrasled;


    for (int j = 1; j <= k; j++) {
        ARototrasled = A;
        BRototrasled = B;
        Eigen::Affine3d* rotoTraslation = rotationAlgorith->rotate(M1, M2);
        calculateRotation(ARototrasled, BRototrasled, rotoTraslation);

        //Calculate the distance between the points

        N.clear();

        for (int i = 0; i < n; i++) {
            distance = (ARototrasled.col(vectorSet[i].first) - BRototrasled.col(vectorSet[i].second)).squaredNorm();
            distance = sqrt(distance);
            threshold = ((j * d) / k);
            if (distance < threshold) {
                N.push_back(vectorSet[i]);
            }

        }


        Eigen::Matrix3Xd N1(3, N.size());
        Eigen::Matrix3Xd N2(3, N.size());
        for (unsigned int i = 0; i < N.size(); i++) {
            N1.col(i) = A.col(N[i].first);
            N2.col(i) = B.col(N[i].second);
        }
        M1 = N1;
        M2 = N2;
    }

    Eigen::Affine3d* rotoTraslation = rotationAlgorith->rotate(M1, M2);
    calculateRotation(A, B, rotoTraslation);

    //Calculate the distance between the points
    M = N;
    for (std::vector<std::pair<int, int> >::iterator it = M.begin(); it != M.end(); ++it) {
        //for (unsigned int i = 0; i < N.size(); i++) {
        distance = (A.col(it->first) - B.col(it->second)).squaredNorm();
        distance = sqrt(distance);
        if (distance > d && E == 'y') {
            M.erase(it);
            if (it == M.end()) {
                break;
            }
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