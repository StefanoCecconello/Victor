#ifndef ROTATOR_H
#define ROTATOR_H

#include <Eigen/Geometry>
#include <Spacer.h>

namespace Victor {
    namespace Biopool {

        class Rotator {
        public:
            // CONSTRUCTORS/DESTRUCTOR:
            Rotator();

            // PREDICATES:
            virtual Eigen::Affine3d* rotate(Eigen::Matrix3Xd set1Matrix, Eigen::Matrix3Xd set2Matrix) = 0;

        

        };

    }
} // namespace

#endif
