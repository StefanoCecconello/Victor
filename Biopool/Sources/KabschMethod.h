#ifndef KABSCHMETHOD_H
#define KABSCHMETHOD_H

#include <Rotator.h>


namespace Victor {
    namespace Biopool {

        class KabschMethod:public Rotator{
        public:

            // CONSTRUCTORS/DESTRUCTOR:
            KabschMethod();
            
            void rotate();

        };

    }
} // namespace

#endif
