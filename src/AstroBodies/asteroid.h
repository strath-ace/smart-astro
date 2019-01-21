#ifndef ASTEROID_H
#define ASTEROID_H

#include "celestial_object.h"

namespace smartastro
{
	namespace astrobodies
	{
		class Asteroid : public Celestial_Object
		{
			public:
				Asteroid();
				
				~Asteroid();
		};
	}
}

#endif // SMARTASTRO_ASTEROID_H
