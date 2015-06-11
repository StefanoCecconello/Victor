#ifndef KABSCHMETHOD_H
#define KABSCHMETHOD_H

#include <Rotator.h>


namespace Victor {
    namespace Biopool {

        class KabschMethod:public Rotator{
        public:

            // CONSTRUCTORS/DESTRUCTOR:
            KabschMethod();
            
            Eigen::Affine3d* rotate(Spacer* set1, Spacer* set2);

        };

    }
} // namespace

#endif
