/*
 * ramd.h
 *
 *  Created on: Apr 10, 2019
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#pragma once

#include <iostream>

#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/ramd_params.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "RandomSphericalDirectionGenerator.h"

struct pull_t;

namespace gmx {

class RAMD
{
public:

	RAMD(RAMDParams const& params, pull_t *pull);

	real add_force(int64_t step, t_mdatoms const& mdatoms,
		gmx::ForceWithVirial *forceWithVirial);

private:

	/// Initialization number for pseudo random number generator
	const RAMDParams params;

	/// Pointer to pull work array
	pull_t *pull;

    /// Random pull direction
	RandomSphericalDirectionGenerator random_spherical_direction_generator;

	/// Current pull direction
	DVec direction;

    /// COM of receptor of last RAMD evaluation step
    DVec com_rec_prev;

    /// COM of ligand of last RAMD evaluation step
    DVec com_lig_prev;

};

} // namespace gmx
