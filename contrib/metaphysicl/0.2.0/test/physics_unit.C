
#include <iostream>

#include "metaphysicl/physics.h"
#include "metaphysicl/numbervector.h"
#include "metaphysicl/numberarray.h"

typedef double Real;

using namespace MetaPhysicL;

DeclareUnaryPhysics(DensityFromSpeciesDensities,
                    DENSITIES_VAR, DENSITY_VAR,
                    sum(rhoi));

DeclareBinaryPhysics(MomentumFromVelocity,
                     DENSITY_VAR, VELOCITY_VAR, MOMENTUM_VAR,
                     U * rho);

DeclareBinaryPhysics(VelocityFromMomentum,
                     DENSITY_VAR, MOMENTUM_VAR, VELOCITY_VAR,
                     rhoU / rho);

DeclareUnaryPhysics(SpeedSquaredFromVelocity,
                    VELOCITY_VAR, SPEED_SQUARED_VAR,
                    U.dot(U));

DeclareBinaryPhysics(SpecificEnergyFromConserved,
                     DENSITY_VAR, ENERGY_VAR, SPECIFIC_ENERGY_VAR,
                     rhoE / rho);

DeclareBinaryPhysics(EnergyFromSpecificEnergy,
                     DENSITY_VAR, SPECIFIC_ENERGY_VAR, ENERGY_VAR,
                     rho * E);

DeclareBinaryPhysics(SpecificInternalEnergyFromSpecificEnergy,
                     SPECIFIC_ENERGY_VAR, SPEED_SQUARED_VAR,
                     SPECIFIC_INTERNAL_ENERGY_VAR,
                     E - UdotU / 2);

DeclareBinaryPhysics(SpecificEnergyFromSpecificInternalEnergy,
                     SPECIFIC_INTERNAL_ENERGY_VAR, SPEED_SQUARED_VAR,
                     SPECIFIC_ENERGY_VAR,
                     e + UdotU / 2);

DeclareBinaryPhysics(LinearTranslationalRotationalEnergyFromTemperature,
                     TEMPERATURE_VAR,
                     TRANSLATIONAL_ROTATIONAL_SPECIFIC_HEAT_VAR,
                     TRANSLATIONAL_ROTATIONAL_ENERGY_VAR,
                     cv_tr * T);

DeclareUnaryPhysics(SpecificInternalEnergyFromOnlyTransationalRotationalEnergy,
                    TRANSLATIONAL_ROTATIONAL_ENERGY_VAR,
                    SPECIFIC_INTERNAL_ENERGY_VAR,
                    e_tr);


template <typename ScalarT,
          typename ScalarC=ScalarT,
          typename VectorV=NumberVector<3,ScalarT>,
          typename VectorRho=NumberVector<13,ScalarT> >
struct TestPhysics
{
  typedef MetaPhysicL::VectorConstructor<
    DensityFromSpeciesDensities,
    VelocityFromMomentum,
    MomentumFromVelocity,
    SpeedSquaredFromVelocity,
    SpecificEnergyFromConserved,
    SpecificInternalEnergyFromSpecificEnergy,
    SpecificEnergyFromSpecificInternalEnergy,
    EnergyFromSpecificEnergy,
    LinearTranslationalRotationalEnergyFromTemperature,
    SpecificInternalEnergyFromOnlyTransationalRotationalEnergy
  >::type AllPhysics;

  typedef typename UIntStructConstructor<
    DENSITIES_VAR, VectorRho,
    VELOCITY_VAR, VectorV,
    TEMPERATURE_VAR, ScalarT
  >::type primitive_vars;

  typedef UIntVectorConstructor<
    DENSITIES_VAR,
    MOMENTUM_VAR,
    ENERGY_VAR
  >::type conserved_vars;

  typedef typename UIntStructConstructor<
    TRANSLATIONAL_ROTATIONAL_SPECIFIC_HEAT_VAR, ScalarC
  >::type constants;

  typedef typename primitive_vars::template Union<constants>::type
    primitive_inputs;

  typedef typename Equations<AllPhysics>::
    SolveList<primitive_inputs,conserved_vars>::type
    primitive_to_conserved;

  typedef typename Equations<AllPhysics>::
    SolveState<primitive_inputs,primitive_to_conserved>::type state;
};


int main(void)
{
  typedef TestPhysics<Real>::state single_state;
  typedef TestPhysics<Real>::primitive_to_conserved single_transformation;

  single_state state1;

  state1.var<DENSITIES_VAR>() = 0; // Don't use uninitialized data!
  state1.var<DENSITIES_VAR>()[0] = 0.1; // species 0 density in kg/m^3
  state1.var<VELOCITY_VAR>() = 0; // Don't use uninitialized data!
  state1.var<VELOCITY_VAR>()[1] = 1000; // y-velocity in m/s
  state1.var<TEMPERATURE_VAR>() = 300; // temp in Kelvin
  state1.var<TRANSLATIONAL_ROTATIONAL_SPECIFIC_HEAT_VAR>() = 750; // c_v in J/kg-K

  single_transformation::ForEach()(EvaluatePhysics<single_state>(state1));

  std::cout << "Densities rho_i = " << state1.var<DENSITIES_VAR>() 
            << " kg/m^3" << std::endl;

  std::cout << "Density rho = " << state1.var<DENSITY_VAR>() 
            << " kg/m^3" << std::endl;

  std::cout << "Speed squared U*U = " << state1.var<SPEED_SQUARED_VAR>() 
            << " m^2/s^2" << std::endl;

  std::cout << "Specific internal energy e = " << state1.var<SPECIFIC_INTERNAL_ENERGY_VAR>() 
            << " J/m^3" << std::endl;

  std::cout << "Specific energy E = " << state1.var<SPECIFIC_ENERGY_VAR>() 
            << " J/m^3" << std::endl;

  std::cout << "Energy rho*E = " << state1.var<ENERGY_VAR>() 
            << " J/m^3" << std::endl;

  typedef NumberArray<3,Real> Scalar2;
  typedef TestPhysics<Scalar2>::state multi_state;
  typedef TestPhysics<Scalar2>::primitive_to_conserved
    multi_transformation;

  multi_state state2;

  state2.var<DENSITIES_VAR>() = 0; // Don't use uninitialized data!
  state2.var<DENSITIES_VAR>()[0][0] = 0.1; // species 0 density in kg/m^3
  state2.var<VELOCITY_VAR>() = 0; // Don't use uninitialized data!
  state2.var<VELOCITY_VAR>()[1][0] = 1000; // y-velocity in m/s
  state2.var<TEMPERATURE_VAR>() = 300; // temp in Kelvin
  state2.var<TRANSLATIONAL_ROTATIONAL_SPECIFIC_HEAT_VAR>() = 750; // c_v in J/kg-K

  multi_transformation::ForEach()(EvaluatePhysics<multi_state>(state2));

  std::cout << "Densities rho_i = " << state2.var<DENSITIES_VAR>() 
            << " kg/m^3" << std::endl;

  std::cout << "Density sums rho = " << sum(state2.var<DENSITIES_VAR>())
            << " kg/m^3" << std::endl;

  std::cout << "Velocities U_i = " << state2.var<VELOCITY_VAR>() 
            << " m/s" << std::endl;

  std::cout << "Temperatures T = " << state2.var<TEMPERATURE_VAR>() 
            << " K" << std::endl;

  std::cout << "Specific Heats c_v = " <<
	       state2.var<TRANSLATIONAL_ROTATIONAL_SPECIFIC_HEAT_VAR>() 
            << " J/kg-K" << std::endl;

  std::cout << "Densities rho = " << state2.var<DENSITY_VAR>() 
            << " kg/m^3" << std::endl;

  std::cout << "Speeds squared U*U = " << state2.var<SPEED_SQUARED_VAR>() 
            << " m^2/s^2" << std::endl;

  std::cout << "Specific internal energies e = " << state2.var<SPECIFIC_INTERNAL_ENERGY_VAR>() 
            << " J/m^3" << std::endl;

  std::cout << "Specific energies E = " << state2.var<SPECIFIC_ENERGY_VAR>() 
            << " J/m^3" << std::endl;

  std::cout << "Energies rho*E = " << state2.var<ENERGY_VAR>() 
            << " J/m^3" << std::endl;
};
