#ifndef KABSCHMETHOD_H
#define KABSCHMETHOD_H

#include <Rotator.h>


namespace Victor {
    namespace Biopool {

        /**@brief Implementation of the kabsch method for find the optimal 
         * rotoslation for superimpose two molecules.
         * 
         * */
        class KabschMethod:public Rotator{
        public:

            // CONSTRUCTORS/DESTRUCTOR:
            KabschMethod();
            
            Eigen::Affine3d* rotate(Eigen::Matrix3Xd set1Matrix, Eigen::Matrix3Xd set2Matrix) const;

        };

    }
} // namespace

#endif
