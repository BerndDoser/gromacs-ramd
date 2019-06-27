/*
 * RandomSphericalDirectionGenerator.h
 *
 *  Created on: Jun 25, 2019
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#pragma once

#include <random>

#include "gromacs/math/vectypes.h"

namespace gmx {

class RandomSphericalDirectionGenerator
{
public:

	RandomSphericalDirectionGenerator(int64_t seed)
     : engine(seed),
       dist(0.0, 1.0)
    {}

	DVec operator()()
	{
	    auto theta = 2 * M_PI * dist(engine);  // azimuthal angle
	    auto psi   =     M_PI * dist(engine);  // polar angle

	    DVec direction;
	    direction[0] = std::cos(theta) * std::sin(psi);
	    direction[1] = std::sin(theta) * std::sin(psi);
	    direction[2] = std::cos(psi);

	    return direction;
	}

private:

	/// Random number generator
	std::default_random_engine engine;

	/// Random number distribution
    std::uniform_real_distribution<> dist;

};

} // namespace gmx
