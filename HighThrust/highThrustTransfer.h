/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATION_PATCHED_CONIC_FULL_PO
#define TUDAT_PROPAGATION_PATCHED_CONIC_FULL_PO

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"

namespace tudat
{

namespace propagators
{

//! Function to setup a body map corresponding to the assumptions of a patched conics trajectory,
//! using default ephemerides for the central and transfer bodies.
/*!
 * Function to setup a body map for the patched conics trajectory. The body map contains the central body, the transfer
 * bodies and the body to be propagated. The positions of the central and transfer bodies are directly retrieved from ephemerides.
 * \param nameCentralBody Name of the central body.
 * \param nameBodyToPropagate Name of the body to be propagated.
 * \param nameTransferBodies Vector containing the names of the transfer bodies.
 * \return Body map for the patched conics trajectory.
 */
simulation_setup::NamedBodyMap setupBodyMapFromEphemeridesForPatchedConicsTrajectoryPO(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& nameTransferBodies );

}

}


#endif // TUDAT_PROPAGATION_PATCHED_CONIC_FULL_PO
