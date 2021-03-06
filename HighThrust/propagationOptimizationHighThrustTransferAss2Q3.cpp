/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <chrono>

#include <Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Astrodynamics/TrajectoryDesign/trajectory.h>
#include <Tudat/SimulationSetup/PropagationSetup/propagationPatchedConicFullProblem.h>
#include <Tudat/SimulationSetup/PropagationSetup/propagationLambertTargeterFullProblem.h>
#include <tudatApplications/AE4866-Assignments-P-O/HighThrust/highThrustTransfer.h>
#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"



#include "../applicationOutput.h"
#include <cmath>

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::gravitation;
using namespace tudat::numerical_integrators;
using namespace tudat::transfer_trajectories;
using namespace tudat::statistics;
using namespace orbital_element_conversions;

//! Function to directly setup a vector of acceleration maps for a patched conics trajectory.
std::vector < basic_astrodynamics::AccelerationMap > getAccelerationModelsPerturbedPatchedConicsTrajectory(
        const double numberOfLegs,
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< std::string >& transferBodyOrder)
{

    std::vector< basic_astrodynamics::AccelerationMap > accelerationMapsVector;
    for (int i = 0 ; i < numberOfLegs ; i++)
    {
        SelectedAccelerationMap accelerationSettingsMap;
        // Acceleration central body

        accelerationSettingsMap[ nameBodyToPropagate ][ nameCentralBody ].push_back(
                    std::make_shared< simulation_setup::AccelerationSettings >(
                        basic_astrodynamics::central_gravity ) );

        // Cannonbal Radiation

        accelerationSettingsMap[ nameBodyToPropagate ][ nameCentralBody ].push_back(
                    std::make_shared< simulation_setup::AccelerationSettings>(
                        basic_astrodynamics::cannon_ball_radiation_pressure));

        // Central Gravity Departure planet
        accelerationSettingsMap[ nameBodyToPropagate ][ transferBodyOrder.at( i ) ].push_back(
                    std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );

        // Central gravity Saturn

        accelerationSettingsMap[ nameBodyToPropagate ][ "Saturn" ].push_back(
                    std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );


        // Central gravity Uranus

        accelerationSettingsMap[ nameBodyToPropagate ][ "Uranus" ].push_back(
                    std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );



        if( i != numberOfLegs -1 )
        {
            // Central gravity Jupiter

            if( transferBodyOrder.at( i ) == "Jupiter" || transferBodyOrder.at(i+1) == "Jupiter" ){
                //if( i < 3){
            }
            else{

                accelerationSettingsMap[ nameBodyToPropagate ][ "Jupiter" ].push_back(
                            std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );

            }
            // Central gravity Mars
            if( transferBodyOrder.at( i ) == "Mars" || transferBodyOrder.at(i+1) == "Mars" ){
                //if( i < 3){
            }
            else{

                accelerationSettingsMap[ nameBodyToPropagate ][ "Mars" ].push_back(
                            std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );

            }
            // Central gravity Earth
            if( transferBodyOrder.at( i ) == "Earth" || transferBodyOrder.at(i+1) == "Earth" ){
                //if( i < 3){
            }
            else{
                accelerationSettingsMap[ nameBodyToPropagate ][ "Earth" ].push_back(
                            std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );

            }
            // Central gravity Venus
            if( transferBodyOrder.at( i ) == "Venus" || transferBodyOrder.at(i+1) == "Venus" ){
                //if( i < 3){
            }
            else{
                accelerationSettingsMap[ nameBodyToPropagate ][ "Venus" ].push_back(
                            std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );

            }

            if( transferBodyOrder.at( i ) != transferBodyOrder.at( i + 1 ) )
            {
                // central gravity Arrival planet
                accelerationSettingsMap[ nameBodyToPropagate ][ transferBodyOrder.at( i + 1 ) ].push_back(
                            std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );

            }

        }

        accelerationMapsVector.push_back( createAccelerationModelsMap(
                                              bodyMap, accelerationSettingsMap, { nameBodyToPropagate }, { nameCentralBody } ) );
    }

    return accelerationMapsVector;

}
/*!
 *   This function computes a patched conic trajectory with a given set of flyby bodies, minimum periapsis distances. The
 *   order of bodies is defined as Earth-Venus-X-Y-Jupiter, with X and Y user-defined. A Trajectory object is created that
 *   computes the required Delta V at arrival and at each flyby (no DSMs are considered). The departure Delta V is not included
 *   in the final Delta V
 *
 *   Subsequently, the function fullPropagationPatchedConicsTrajectory is called, which propagates the problem
 *   numerically in the following way:
 *
 *   - The patched conic (consisting of a set of Lambert targeters), is reconstructed, and the output (vehicle state) at regular
 *     intervals is saved.
 *   - Starting at the midpoint (in time) of each transfer leg, the spacecraft state is numerically propagated using acceleration/
 *     environment settings that may be defined by the user. The propagation is done from the midpoint forwards, and backwards,
 *     to obtain the state history of the full leg.
 *
 *   In the code, as given to you, the dynamical model used in the numerical propagation is not fully consistent with that of the
 *   patched conic method. In addition to the point mass gravity attraction by the Sun (Sun fixed at the origin), the perturbations
 *   of the arrival and departure planet are also taken into account (as point-mass gravities)
 *
 *   Key outputs (per leg):
 *
 *   lambertTargeterResultForEachLeg: a list of the state history of the spacecraft (per leg) according to the patched conic
 *      method
 *   fullProblemResultForEachLeg: a list of the state history of the spacecraft (per leg) as produced by the numerical propagation
 *
 *   Input parameters:
 *
 *   trajectoryIndependentVariables: A vector defining the start time (first entry, in seconds since J2000) and duration of all 4
 *      legs (second, third, fourth and fifth entries respectively). The final entry
 *      of this vector is an integer defining which entry of the transferCaseNames vector is to be used,
 *      e.g. which bodies X and Y are used for the 3rd and 4th body.
 */
int main( )
{

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );


    std::string outputPath = tudat_applications::getOutputPath( "HighThrustAss2Q3" );
    //std::map< double, Eigen::VectorXd >  refMap;
    std::map< int, std::map< double, Eigen::VectorXd > > refMap;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        TRANSFER SETTINGS                 //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define different cases for 3rd and 4th body
    std::vector< std::string > transferCases = { "EVEEJ", "EVVEJ",  "EVEVJ", "EVVMJ", "EVEMJ", "EVMMJ", "EVMVJ"  };
    std::vector< std::pair< std::string, std::string > > transferCaseNames =
    { { "Earth", "Earth" }, { "Venus", "Earth" }, { "Earth", "Venus" }, { "Venus", "Mars" }, { "Earth", "Mars" },
      { "Mars", "Mars" }, { "Mars", "Venus" } };

    // DEFINE PROBLEM INDEPENDENT VARIABLES HERE:
    std::vector< double > trajectoryParameters =
    { -5185.040933880841,	 284.508757767517, 	 193.7914543558609,	 282.9403211987581,	 767.7274155863713,	6 };

    int transferCase = trajectoryParameters.at( 5 );

    // Set body order (no DSM) for current settings
    std::vector< std::string > transferBodyOrder =
    { "Earth", "Venus", transferCaseNames.at( transferCase ).first, transferCaseNames.at( transferCase ).second, "Jupiter" };
    std::vector< TransferLegType > transferLegTypes = { mga_Departure, mga_Swingby, mga_Swingby, mga_Swingby, capture };


    // Define settings for capture at target planet
    double captureSemiMajorAxis = 1.0895e8 / 0.02;
    double captureEccentricity = 0.98;
    std::vector< double > departureCaptureSemiMajorAxes = { TUDAT_NAN, captureSemiMajorAxis };
    std::vector< double > departureCaptureEccentricities = { TUDAT_NAN, captureEccentricity };
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        SETUP SOLAR SYSTEM BODIES            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create body map
    NamedBodyMap bodyMapForPatchedConic = setupBodyMapFromEphemeridesForPatchedConicsTrajectoryPO(
                "Sun", "Spacecraft", transferBodyOrder );

    bodyMapForPatchedConic[ "Saturn" ] = std::make_shared< Body >( );
    bodyMapForPatchedConic[ "Saturn" ]->setEphemeris(
                std::make_shared< ephemerides::ApproximatePlanetPositions>(
                    "Saturn") );
    bodyMapForPatchedConic[ "Saturn" ]->setGravityFieldModel(
                createGravityFieldModel(
                    std::make_shared< CentralGravityFieldSettings >(
                        spice_interface::getBodyGravitationalParameter(
                            "Saturn"  ) ),  "Saturn"  ) );

    bodyMapForPatchedConic[ "Uranus" ] = std::make_shared< Body >( );
    bodyMapForPatchedConic[ "Uranus" ]->setEphemeris(
                std::make_shared< ephemerides::ApproximatePlanetPositions>(
                    "Uranus") );
    bodyMapForPatchedConic[ "Uranus" ]->setGravityFieldModel(
                createGravityFieldModel(
                    std::make_shared< CentralGravityFieldSettings >(
                        spice_interface::getBodyGravitationalParameter(
                            "Uranus"  ) ),  "Uranus"  ) );








    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE SPACECRAFT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMapForPatchedConic[ "Spacecraft" ] = std::make_shared< simulation_setup::Body >( );
    bodyMapForPatchedConic[ "Spacecraft" ]->setConstantBodyMass( 4362.0);

    double referenceAreaRadiation = 51; // AANPASSEN
    double radiationPressureCoefficient = 1.2; // AANPASSEN
    std::shared_ptr< RadiationPressureInterfaceSettings > JUPRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient);

    // Create and set radiation pressure settings
    bodyMapForPatchedConic[ "Spacecraft" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    JUPRadiationPressureSettings, "Spacecraft", bodyMapForPatchedConic ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMapForPatchedConic, "SSB", "ECLIPJ2000" );


    unsigned int prange;
    unsigned int erange;
    std::vector<double> errorList = {1000,1000,1000,0.01,0.01,0.01};
    for( unsigned int k = 0; k <7; k++){
        std::map< int ,Eigen::Vector4d> MaxPosError;
        std::map< int ,Eigen::Vector4d> MaxVelError;
        int x = 0;

        if(k == 0){
            prange = 1;
            erange = 1;
        }
        else {
            prange = 1000 ; // prob 500
            erange = 5;
        }
        for( unsigned int e = 0; e<erange; e++){
            double seed = 0.0;
            std::function<double()> randomNumberGenerator;
            if(k > 0){
                randomNumberGenerator = createBoostContinuousRandomVariableGeneratorFunction( normal_boost_distribution, {0.0, errorList.at(k-1)},seed);
            }
            seed++;
            for( unsigned int p = 0; p<prange; p++){
                double error;




                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////             CREATE PATCHED CONIC SEMI-ANALYTICAL TRAJECTORY            ///////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                // Get list of minimum flyby periapsis radii
                std::vector< double > minimumPericenterRadii = getDefaultMinimumPericenterRadii(
                            transferBodyOrder );

                std::vector< double > trajectoryIndependentVariables;
                for( unsigned int i = 0; i < trajectoryParameters.size( ) - 1; i++ )
                {        trajectoryIndependentVariables.push_back( trajectoryParameters.at( i ) * physical_constants::JULIAN_DAY );
                }

                // Add entry for interface consistency
                trajectoryIndependentVariables.push_back( TUDAT_NAN );

                // Create patched conic calculation object (no numerical propagation; departure Delta V not included)
                transfer_trajectories::Trajectory trajectory = createTransferTrajectoryObject(
                            bodyMapForPatchedConic, transferBodyOrder, "Sun", transferLegTypes, trajectoryIndependentVariables,
                            minimumPericenterRadii, false, TUDAT_NAN, TUDAT_NAN, true,
                            captureSemiMajorAxis, captureEccentricity );


                // Retrieve total Delta V values
                double totalDeltaV;
                trajectory.calculateTrajectory( totalDeltaV );
                double captureDeltaV;
                trajectory.getCaptureDeltaV( captureDeltaV );

                std::cout<<"Total/capture Delta V: "<<totalDeltaV<<" "<<captureDeltaV<<std::endl;

                // Retrieve times, positions, and delta V at each maneuver
                std::vector < Eigen::Vector3d > positionVector;
                std::vector < double > timeVector;
                std::vector < double > deltaVVector;
                trajectory.maneuvers( positionVector, timeVector, deltaVVector );


                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////             NUMERICALLY PROPAGATE DYNAMICS            ////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                // Define environment for propagation (equal to that of patched conic)
                NamedBodyMap bodyMapForPropagation = bodyMapForPatchedConic;

                // Define acceleration settings
                std::vector< basic_astrodynamics::AccelerationMap > accelerationMap =
                        getAccelerationModelsPerturbedPatchedConicsTrajectory(
                            transferLegTypes.size( ), "Sun", "Spacecraft", bodyMapForPropagation, transferBodyOrder);

                // Define integrator settings
                double relativeTolerance = 1E-9;
                double absoluteTolerance = 1E-9;
                double minimumStepSize   = std::numeric_limits< double >::epsilon( );
                double maximumStepSize   = std::numeric_limits< double >::infinity( );
                double initialStepSize   = 1000;
                double initialTime = TUDAT_NAN;
                std::shared_ptr< IntegratorSettings < double > > integratorSettings =  std::make_shared < RungeKuttaVariableStepSizeSettings <double> > (
                            initialTime, initialStepSize, RungeKuttaCoefficients::rungeKuttaFehlberg78, minimumStepSize,
                            maximumStepSize, relativeTolerance, absoluteTolerance );

                // Create list of relevant bodies
//                std::vector< std::string > bodyList;
//                for( unsigned int i = 0; i < transferBodyOrder.size( ); i++ )
//                {
//                    if( std::find( bodyList.begin( ), bodyList.end( ), transferBodyOrder.at( i ) ) ==
//                            bodyList.end( ) )
//                    {
//                        bodyList.push_back( transferBodyOrder.at( i ) );
//                    }
//                }
//                bodyList.push_back( "Sun" );

//                // Create list of dependent variables to save (distance to all flyby bodies and Sun)
//                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariableList;
//                for( unsigned int i = 0; i < bodyList.size( ); i++ )
//                {
//                    dependentVariableList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
//                                                         relative_distance_dependent_variable, "Spacecraft", bodyList.at( i ) ) );

//                }
                // Save dependent variables for each propagation leg
                std::vector< std::shared_ptr< DependentVariableSaveSettings > > dependentVariablesToSave;
//                for( unsigned int j = 0; j < transferBodyOrder.size( ); j++ )
//                {
//                    dependentVariablesToSave.push_back( std::make_shared< DependentVariableSaveSettings >(
//                                                            dependentVariableList ) );
//                }

                // Define propagator type
                TranslationalPropagatorType propagatorType = cowell;

                // Create propagator settings for patched conic (per arc; backward and forward from arc midpoint)
                // Propagation currently terminates on sphere of influence of body.
                std::vector< std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
                        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > > propagatorSettings =
                        getPatchedConicPropagatorSettings(
                            bodyMapForPropagation, accelerationMap, transferBodyOrder, "Sun", "Spacecraft", transferLegTypes,
                            trajectoryIndependentVariables, minimumPericenterRadii, departureCaptureSemiMajorAxes,
                            departureCaptureEccentricities, dependentVariablesToSave, propagatorType, true );

                // Create an object of `steady_clock` class
                std::chrono::steady_clock sc;

                // Start timer
                auto start = sc.now();

                // Propagate full dynamics of problem
                std::map< int, std::map< double, Eigen::Vector6d > > lambertTargeterResultForEachLeg;
                std::map< int, std::map< double, Eigen::Vector6d > > fullProblemResultForEachLeg;
                std::map< int, std::map< double, Eigen::VectorXd > > dependentVariableResultForEachLeg;
                fullPropagationPatchedConicsTrajectory(
                            bodyMapForPropagation, transferBodyOrder,
                            "Sun", transferLegTypes, trajectoryIndependentVariables, minimumPericenterRadii,
                            departureCaptureSemiMajorAxes, departureCaptureEccentricities,
                            propagatorSettings, integratorSettings,
                            lambertTargeterResultForEachLeg, fullProblemResultForEachLeg, dependentVariableResultForEachLeg );

                // End timer (starting & ending is done by measuring the time at the moment the process started & ended respectively)
                auto end = sc.now();


                // Measure time span between start & end
                auto time_span = static_cast<std::chrono::duration<double>>(end - start);
                double runTimeInSeconds = time_span.count( );
                //result.insert(std::pair<unsigned int,double>(k,runTimeInSeconds));
                //                // Write perturbed satellite propagation history to file.
                //                input_output::writeDataMapToTextFile( result, "timeResult"+std::to_string(k)+std::to_string(e)+std::to_string(p) + ".dat", outputPath );

                std::cout<<"Operation took: "<<runTimeInSeconds<<" seconds"<<std::endl;


                std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
                        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 );
                std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > > benchmarkInterpolator;
                std::map < double, Eigen::VectorXd > interpolatedState;

                std::vector< double > positionError;
                std::vector< double > velocityError;
                Eigen::Vector4d posErr(0,0,0,0);
                Eigen::Vector4d velErr(0,0,0,0);
                if ( k == 0){

                    for( auto resultIterator : fullProblemResultForEachLeg )
                    {

                        int currentArc = resultIterator.first;
                        double currentArcMiddleTime = trajectoryParameters.at( 0 ) + trajectoryParameters.at( 1 ) / 2.0;
                        // Retrieve state history for current arc
                        std::map< double, Eigen::Vector6d > fullProblemSolution = resultIterator.second;

                        // Retrieve numerical state at middle of arc.
                        Eigen::Vector6d currentArcMiddleState = interpolators::createOneDimensionalInterpolator(
                                    fullProblemSolution, std::make_shared< interpolators::LagrangeInterpolatorSettings >(
                                        8 ) )->interpolate( currentArcMiddleTime * 86400.0 );

                        // Reset integrator initial time
                        integratorSettings->initialTime_ = currentArcMiddleTime * 86400.0;

                        // Retrieve propagation settings for forward propagation, and reset initial state/final time
                        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > forwardPropagatorSettings =
                                propagatorSettings.at( resultIterator.first ).second;
                        forwardPropagatorSettings->resetInitialStates( currentArcMiddleState );
                        forwardPropagatorSettings->resetTerminationSettings( std::make_shared< PropagationTimeTerminationSettings >(
                                                                                 fullProblemSolution.rbegin( )->first ) );

                        // Ensure time step is positive (forward integration)
                        integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

                        // Propagate dynamics forward and print results to file
                        SingleArcDynamicsSimulator< > forwardDynamicsSimulator(
                                    bodyMapForPropagation, integratorSettings, forwardPropagatorSettings );
                        std::map< double, Eigen::VectorXd > integrationResult;
                        integrationResult = forwardDynamicsSimulator.getEquationsOfMotionNumericalSolution( );




                        // Update arc middle time for next arc.

                        if( currentArc < fullProblemResultForEachLeg.size( ) - 1 )
                        {
                            currentArcMiddleTime += ( trajectoryParameters.at( currentArc + 1 ) + trajectoryParameters.at( currentArc + 2 ) ) / 2.0;

                        }
                        refMap.insert(std::pair<int, std::map< double, Eigen::VectorXd >>(resultIterator.first,integrationResult));
                        //refMap=integrationResult;
                    }
                }
                if ( k != 0){
                    int i = 0;
                    for( auto resultIterator : fullProblemResultForEachLeg )
                    {
                        std::vector< double> PositionReference;
                        std::vector< double> VelocityReference;
                        std::vector< double> PositionRun;
                        std::vector< double> VelocityRun;
                        std::map < double, Eigen::Vector6d > interpolatedState;
                        std::map< double, Eigen::Vector6d> resulPerLag = resultIterator.second;

                        int currentArc = resultIterator.first;
                        double currentArcMiddleTime = trajectoryParameters.at( 0 ) + trajectoryParameters.at( 1 ) / 2.0;
                        // Retrieve state history for current arc
                        std::map< double, Eigen::Vector6d > fullProblemSolution = resultIterator.second;

                        // Retrieve numerical state at middle of arc.
                        Eigen::Vector6d currentArcMiddleState = interpolators::createOneDimensionalInterpolator(
                                    fullProblemSolution, std::make_shared< interpolators::LagrangeInterpolatorSettings >(
                                        8 ) )->interpolate( currentArcMiddleTime * 86400.0 );

                        // Reset integrator initial time
                        integratorSettings->initialTime_ = currentArcMiddleTime * 86400.0;
                        currentArcMiddleState[k-1] = currentArcMiddleState[k-1] + randomNumberGenerator();

                        // Retrieve propagation settings for forward propagation, and reset initial state/final time
                        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > forwardPropagatorSettings =
                                propagatorSettings.at( resultIterator.first ).second;
                        forwardPropagatorSettings->resetInitialStates( currentArcMiddleState );
                        forwardPropagatorSettings->resetTerminationSettings( std::make_shared< PropagationTimeTerminationSettings >(
                                                                                 fullProblemSolution.rbegin( )->first ) );

                        // Ensure time step is positive (forward integration)
                        integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );


                        // Propagate dynamics forward and print results to file
                        SingleArcDynamicsSimulator< > forwardDynamicsSimulator(
                                    bodyMapForPropagation, integratorSettings, forwardPropagatorSettings );
                        std::map< double, Eigen::VectorXd > integrationResult;
                        integrationResult = forwardDynamicsSimulator.getEquationsOfMotionNumericalSolution( );

                        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > benchmarkInterpolator =
                                interpolators::createOneDimensionalInterpolator(
                                    refMap.at(i), interpolatorSettings );
                        for( auto timeIterator : integrationResult){
                            interpolatedState[ timeIterator.first ] =
                                    benchmarkInterpolator->interpolate(
                                        timeIterator.first );
                        }
                        i++;

                        for( auto timeIterator : interpolatedState){
                            PositionReference.push_back(sqrt(pow(timeIterator.second(0),2.0)+pow(timeIterator.second(1),2.0)+pow(timeIterator.second(2),2.0)));
                            VelocityReference.push_back(sqrt(pow(timeIterator.second(3),2.0)+pow(timeIterator.second(4),2.0)+pow(timeIterator.second(5),2.0)));
                        }
                        for( auto timeIterator : integrationResult){
                            PositionRun.push_back(sqrt(pow(timeIterator.second(0),2.0)+pow(timeIterator.second(1),2.0)+pow(timeIterator.second(2),2.0)));
                            VelocityRun.push_back(sqrt(pow(timeIterator.second(3),2.0)+pow(timeIterator.second(4),2.0)+pow(timeIterator.second(5),2.0)));
                        }
                        std::vector< double> differencePos;
                        for( unsigned int w = 4; w < PositionReference.size()-4; w++){
                            double pos = PositionReference.at(w)-PositionRun.at(w);
                            differencePos.push_back((pos));
                        }
                        positionError.push_back(*max_element(differencePos.begin(), differencePos.end()));

                        std::vector< double> differenceVel;
                        for( unsigned int w = 4; w <VelocityReference.size()-4; w++){
                            double vel = VelocityReference.at(w)-VelocityRun.at(w);
                            differenceVel.push_back((vel));
                        }
                        velocityError.push_back( *max_element(differenceVel.begin(), differenceVel.end()));
                        std::cout<<"------------------------"<< velocityError.size() << " "<<positionError.size()<<std::endl;;


                        // Update arc middle time for next arc.

                        if( currentArc < fullProblemResultForEachLeg.size( ) - 1 )
                        {
                            currentArcMiddleTime += ( trajectoryParameters.at( currentArc + 1 ) + trajectoryParameters.at( currentArc + 2 ) ) / 2.0;

                        }


                    }
                    for( unsigned a = 0; a<4; a++){
                        posErr(a) = positionError.at(a);
                        velErr(a) = velocityError.at(a);
                    }
                }
                MaxPosError.insert( std::pair<int,Eigen::Vector4d >(x , posErr));
                MaxVelError.insert( std::pair<int,Eigen::Vector4d >(x , velErr));
                x++;
                //                    //        if( k == 2){
                //                    //        double bodyMass = bodyMapForPropagation[ "Spacecraft"]->getBodyMass();
                //                    //        double deltaV = deltaVVector.at(currentArc);
                //                    //        double Isp = 300;
                //                    //        double g0 = 9.80665;
                //                    //        double newBodyMass = bodyMass*exp(-deltaV/Isp/g0);
                //                    //        bodyMapForPropagation[ "Spacecraft" ]->setConstantBodyMass(newBodyMass);
                //                    //        //bodyMapForPropagation[ "Spacecraft" ]->updateMass();
                //                    //        double mass = bodyMapForPropagation[ "Spacecraft"]->getBodyMass();
                //                    //        std::cout<<"mass: "<<mass<<std::endl;
                //                    //        }

                //                }

                //                // Write patched conic results to file for each leg
                //                for( auto resultIterator : lambertTargeterResultForEachLeg )
                //                {
                //                    input_output::writeDataMapToTextFile(
                //                                resultIterator.second, "lambertResult" +
                //                                std::to_string( resultIterator.first ) +std::to_string(k)+std::to_string(e)+std::to_string(p)+ ".dat", outputPath );
                //                }

                //                // Write numerical propagation results to file for each leg
                //                for( auto resultIterator : fullProblemResultForEachLeg )
                //                {

                //                    input_output::writeDataMapToTextFile(
                //                                resultIterator.second, "numericalResult" +
                //                                std::to_string( resultIterator.first )+std::to_string(k)+std::to_string(e)+std::to_string(p) + ".dat", outputPath );
                //                }

                //                // Write numerical propagation results to file for each leg
                //                for( auto resultIterator : dependentVariableResultForEachLeg )
                //                {

                //                    input_output::writeDataMapToTextFile(
                //                                resultIterator.second, "dependentResult" +
                //                                std::to_string( resultIterator.first )+ std::to_string(k)+std::to_string(e)+std::to_string(p) + ".dat", outputPath );

                //                }
                //spice_interface::clearSpiceKernels();
            }

        }
        input_output::writeDataMapToTextFile(
                    MaxPosError, "PositionResultQ3PerScenario" +
                    std::to_string(k)+".dat", outputPath );
        input_output::writeDataMapToTextFile(
                    MaxVelError, "VelocityResultQ3PerScenario" +
                    std::to_string(k)+".dat", outputPath );
    }

    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
