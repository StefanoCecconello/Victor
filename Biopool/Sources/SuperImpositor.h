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
            double calculateMaxSub();



            Spacer* getSet2() const;
            Spacer* getSet1() const;
            Spacer getRMSDset2() const;
            Spacer getRMSDset1() const;

        private:

            // PREDICATES:
            void calculateRotation(Eigen::Matrix3Xd& firstSet, Eigen::Matrix3Xd& secondSet, Eigen::Affine3d* rotoTraslation);

            void rotationApplication(Eigen::Matrix3Xd R);
            void translationApplication(Eigen::Vector3d S);

            //Help function
            vgMatrix3<double> fromMatrix3XdTovgMatrix3(Eigen::Matrix3Xd matrix3Xd) const;
            Spacer fromMatrix3XdToSpacer(Eigen::Matrix3Xd matrix3Xd, int num) const;
            vgVector3<double> fromVector3dTovgVector3(Eigen::Vector3d Vector3d) const;
            Eigen::Matrix3Xd fromSpacerToMatrix3Xd(Spacer spacerSet) const;
            std::vector<std::pair<int, int> > maxSubAlignment(Eigen::Matrix3Xd& firstSet, Eigen::Matrix3Xd& secondSet, std::vector< std::pair<int, int> > vectorSet, double d);
            std::vector< std::pair<int, int> > Extend(std::vector<std::pair<int, int> > M, std::vector<std::pair<int, int> > vectorSet, Eigen::Matrix3Xd& A, Eigen::Matrix3Xd& B, double d, int L,int n);


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

