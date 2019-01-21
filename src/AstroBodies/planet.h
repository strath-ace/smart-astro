#ifndef PLANET_H
#define PLANET_H

extern "C" {
#include "../../../CSpice/cspice/include/SpiceUsr.h"
#include "../../../CSpice/cspice/include/SpiceZfc.h"
}

#include "celestial_object.h"

namespace smartastro
{
	namespace astrobodies
	{

		class Planet : public Celestial_Object
		{
			public:
				Planet();

				~Planet();
		};
	}
}

#endif // SMARTASTRO_PLANET_H
