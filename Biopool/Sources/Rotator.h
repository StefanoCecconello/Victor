#ifndef ROTATOR_H
#define ROTATOR_H

#include <Eigen/Geometry>


namespace Victor {
    namespace Biopool {

        class Rotator {
        public:
            // CONSTRUCTORS/DESTRUCTOR:
            Rotator();
            
            // PREDICATES:
            virtual void rotate()=0;
            
        private:
            // ATTRIBUTES:
            int set1;
            int set2;
            
        };

    }
} // namespace

#endif
