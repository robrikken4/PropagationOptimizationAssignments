
/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/propagationPatchedConicFullProblem.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationLambertTargeterFullProblem.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"


#include <Tudat/InputOutput/basicInputOutput.h>
#include <tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/applicationOutput.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>


#include "Tudat/Astrodynamics/TrajectoryDesign/captureLeg.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/departureLegMga.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/departureLegMga1DsmPosition.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/departureLegMga1DsmVelocity.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/planetTrajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/swingbyLegMga.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/swingbyLegMga1DsmPosition.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/swingbyLegMga1DsmVelocity.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"
#include "Tudat/Mathematics/BasicMathematics/legendrePolynomials.h"

namespace tudat
{

namespace propagators
{
using namespace ephemerides;
using namespace interpolators;
using namespace numerical_integrators;
using namespace spice_interface;
using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace orbital_element_conversions;
using namespace propagators;
using namespace aerodynamics;
using namespace basic_mathematics;
using namespace input_output;
using namespace estimatable_parameters;


//! Function to setup a body map corresponding to the assumptions of a patched conics trajectory,
//! using default ephemerides for the central and transfer bodies.
simulation_setup::NamedBodyMap setupBodyMapFromEphemeridesForPatchedConicsTrajectoryPO(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& nameTransferBodies)
{
    spice_interface::loadStandardSpiceKernels( );

    // Create central and transfer bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( nameCentralBody );
    for ( unsigned int i = 0 ; i < nameTransferBodies.size( ) ; i ++ )
    {
        bodiesToCreate.push_back( nameTransferBodies[ i ] );
    }


    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ nameCentralBody ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ nameCentralBody ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ nameCentralBody ]->rotationModelSettings->resetOriginalFrame( frameOrientation );

    double sunNormalizedJ2 = 2.0E-7 / calculateLegendreGeodesyNormalizationFactor( 2, 0 );
    double gravitationalParameter = spice_interface::getBodyGravitationalParameter("Sun");
    double averageRadius = spice_interface::getAverageRadius("Sun");

    bodySettings[ "Sun" ]->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                gravitationalParameter, averageRadius,
                ( Eigen::Matrix3d( ) << 1.0, 0.0, 0.0,
                  0.0, 0.0, 0.0,
                  sunNormalizedJ2 , 0.0, 0.0 ).finished( ),
                Eigen::Matrix3d::Zero( ), "IAU_Sun" );

    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );


    // Define body to propagate.
    bodyMap[ nameBodyToPropagate ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ nameBodyToPropagate ]->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                      std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                      < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );


    return bodyMap;
}


}

}

