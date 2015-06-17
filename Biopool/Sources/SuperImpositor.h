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


            Spacer* getSet2() const;
            Spacer* getSet1() const;

        private:
            
            // PREDICATES:
            void rotationApplication(Eigen::Matrix3Xd R);
            void translationApplication(Eigen::Vector3d S);
            Eigen::Matrix3Xd fromSpacerToMatrix3Xd(Spacer spacerSet) const;
            vgMatrix3<double> fromMatrix3XdTovgMatrix3(Eigen::Matrix3Xd matrix3Xd) const;
            vgVector3<double> fromVector3dTovgVector3(Eigen::Vector3d Vector3d) const;
            
            
            // ATTRIBUTES:
            // This is the rotation algorithm chose for this superImpositor
            Rotator* rotationAlgorith;
            //This is the Spacer of the first protein insered
            Spacer* set1;
            //This is the Spacer of the second protein insered
            Spacer* set2;
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

