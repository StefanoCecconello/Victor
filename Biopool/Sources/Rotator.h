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
            virtual double rotate(Spacer* set1, Spacer* set2) = 0;

        

        };

    }
} // namespace

#endif
