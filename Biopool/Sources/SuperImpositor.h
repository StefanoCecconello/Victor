/* 
 * File:   superImpositor.h
 * Author: cecco
 *
 * Created on June 8, 2015, 4:22 PM
 */

#ifndef SUPERIMPOSITOR_H
#define	SUPERIMPOSITOR_H

#include "Protein.h"
#include <KabschMethod.h>
#include "Eigen/Geometry"
#include <vector>

namespace Victor {
    namespace Biopool {

        class SuperImpositor {
        public:
            /**
             * @brief Do the superimposition between two proteins.
             * 
             *  
             *    Ulteriore spiegazione
             * */

            // CONSTRUCTORS/DESTRUCTOR:
            SuperImpositor(Protein* firstProtein, Protein* secondProtein, string method);
            virtual ~SuperImpositor();

            // PREDICATES:
            double calculateRMSD();
            double calculateMaxSub(double d, std::vector<std::pair<int, int> > vectorSet, char E);
            double calculateGdt(std::vector<std::pair<int, int> > vectorSet);
            double calculateTMScore(std::vector<std::pair<int, int> > vectorSet);



            Spacer* getSet2() const;
            Spacer* getSet1() const;
            Spacer getRMSDset2() const;
            Spacer getRMSDset1() const;
            Spacer getMaxSubset2() const;
            Spacer getMaxSubset1() const;
            Spacer getTMScoreset2() const;
            Spacer getTMScoreset1() const;
            Spacer getGdtset4_2() const;
            Spacer getGdtset4_1() const;
            Spacer getGdtset3_2() const;
            Spacer getGdtset3_1() const;
            Spacer getGdtset2_2() const;
            Spacer getGdtset2_1() const;
            Spacer getGdtset1_2() const;
            Spacer getGdtset1_1() const;

        private:

            // PREDICATES:
            double maxEvaluate(double d, std::vector<std::pair<int, int> > vectorSet, char E, Eigen::Matrix3Xd& modifyMatrixSet1, Eigen::Matrix3Xd& modifyMatrixSet2, Eigen::Affine3d* rotoTraslation);
            void calculateRotation(Eigen::Matrix3Xd& firstSet, Eigen::Matrix3Xd& secondSet, Eigen::Affine3d* rotoTraslation);

            void rotationApplication(Eigen::Matrix3Xd R);
            void translationApplication(Eigen::Vector3d S);

            //Help function
            vgMatrix3<double> fromMatrix3XdTovgMatrix3(Eigen::Matrix3Xd matrix3Xd) const;
            Spacer fromMatrix3XdToSpacer(Eigen::Matrix3Xd matrix3Xd, int num) const;
            Spacer rotateSpacer(Eigen::Affine3d* rotoTraslation,int num) const;
            vgVector3<double> fromVector3dTovgVector3(Eigen::Vector3d Vector3d) const;
            Eigen::Vector3d fromvgVector3ToVector3d(vgVector3<double> vgVector3) const;
            Eigen::Matrix3Xd fromSpacerToMatrix3Xd(Spacer spacerSet) const;
            std::vector<std::pair<int, int> > maxSubAlignment(Eigen::Matrix3Xd& firstSet, Eigen::Matrix3Xd& secondSet, std::vector< std::pair<int, int> > vectorSet, double d, char E, Eigen::Affine3d* rotoTraslation);
            std::vector< std::pair<int, int> > Extend(std::vector<std::pair<int, int> > M, std::vector<std::pair<int, int> > vectorSet, Eigen::Matrix3Xd& A, Eigen::Matrix3Xd& B, double d, int L, int n, char E, Eigen::Affine3d* rotoTraslation);


            // ATTRIBUTES:
            // This is the rotation algorithm chose for this superImpositor
            Rotator* rotationAlgorith;
            //This is the Spacer of the first protein insered
            Spacer* set1;
            //This is the Spacer of the second protein insered
            Spacer* set2;
            //This is the Spacer of the first protein insered for  RMSD method
            Spacer RMSDset1;
            //This is the Spacer of the second protein insered for  RMSD method
            Spacer RMSDset2;
            //This is the Spacer of the first protein insered for  RMSD method
            Spacer MaxSubset1;
            //This is the Spacer of the second protein insered for  RMSD method
            Spacer MaxSubset2;
            
            //This is the Spacer of the first protein insered for  RMSD method
            Spacer Gdtset1_1;
            //This is the Spacer of the second protein insered for  RMSD method
            Spacer Gdtset1_2;
            //This is the Spacer of the first protein insered for  RMSD method
            Spacer Gdtset2_1;
            //This is the Spacer of the second protein insered for  RMSD method
            Spacer Gdtset2_2;
            //This is the Spacer of the first protein insered for  RMSD method
            Spacer Gdtset3_1;
            //This is the Spacer of the second protein insered for  RMSD method
            Spacer Gdtset3_2;
            //This is the Spacer of the first protein insered for  RMSD method
            Spacer Gdtset4_1;
            //This is the Spacer of the second protein insered for  RMSD method
            Spacer Gdtset4_2;
            
            
            //This is the Spacer of the first protein insered for  RMSD method
            Spacer TMScoreset1;
            //This is the Spacer of the second protein insered for  RMSD method
            Spacer TMScoreset2;
            //This is a matrix that contain all the coords of the CA atoms in set1
            Eigen::Matrix3Xd matrixSet1;
            //This is a matrix that contain all the coords of the CA atoms in set2
            Eigen::Matrix3Xd matrixSet2;


            //This is the value of the RMSD between the sequences set1 and set2
            double rmsdValue;
            //This is the value of the maxsub between the sequences set1 and set2
            double maxsubValue;
            //This is the value of the gdt between the sequences set1 and set2
            double gdtValue;
            //This is an array that contein the alignment between the sequences set1 and set2 using rmsd method
            int* rmsdAlignment;
            //This is an array that contein the alignment between the sequences set1 and set2 using maxsub method
            int* maxsubAlignment;
            //This is an array that contein the alignment between the sequences set1 and set2 using gdt method
            int* gdtAlignment;
        };
    }
}

#endif	/* SUPERIMPOSITOR_H */

