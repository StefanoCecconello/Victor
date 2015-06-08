#ifndef KABSCHMETHOD_H
#define KABSCHMETHOD_H

#include <Rotator.h>


namespace Victor {
    namespace Biopool {

        class KabschMethod:public Rotator{
        public:

            // CONSTRUCTORS:

            /// Constructor.
            KabschMethod();
            
            void rotate();

        };

    }
} // namespace

#endif
