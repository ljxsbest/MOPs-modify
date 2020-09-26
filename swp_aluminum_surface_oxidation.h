/***************************************************
this file is to solve surface oxidation before gas-phase combustion.

the basic mechanism: molar of O2 consume = alpha*v*dt*S.
where v is calculated by Boltzmann equation.

****************************************************/
#ifndef SWP_ALUMINUM_SURFACE_OXIDATION_H_
#define SWP_ALUMINUM_SURFACE_OXIDATION_H_

#include "swp_surface_reaction.h"
#include "swp_primary.h"

namespace Sweep{
namespace Processes{

class AlSurfOxidation: public Sweep::Processes::SurfaceReaction
{
public:
	// Mechanism constructor
	AlSurfOxidation(const Sweep::Mechanism &mech);

	// Copy constructor
	AlSurfOxidation(const AlSurfOxidation &copy);

	//stream reading constructor
	AlSurfOxidation(std::istream &in,
	const Sweep::Mechanism &mech);

	// Create a copy of the particle process.
	AlSurfOxidation *const AlSurfOxidation::Clone(void) const;

	// Assignment operator
	AlSurfOxidation &operator=(const AlSurfOxidation &rhs);

	// return the process type.
	ProcessType ID(void) const;

	double SurfaceRateConstant(double t,
							const Cell &sys,
							const Particle &sp) const;
	double O2ConsumeRate(double t,
							const Cell &sys,
							const Particle &sp) const;
	double AlConsumeRate(double t,
							const Cell &sys,
							const Particle &sp) const;
	double Al2O3ProduceRate(double t,
							 const Cell &sys,
							 const Particle &sp) const;

	double Rate(double t, const Cell &sys, const Particle &sp) const;

	double HeatProdRate(double t,
					const Cell &sys,
					const Particle &sp) const;

	//writes the object to a binary stream.
	void Serialize(std::ostream &out) const;

	//read the object from a binary stream.
	void Deserialize(
		std::istream &in,
		const Sweep::Mechanism &mech
		);
	double k1;
	double mole_O2Consume;
	double mole_AlConsume;
	double mole_Al2O3Produce;
	double rate_heatprod;

private:
	// Private default constructor
	AlSurfOxidation();

	// Initialize indices for gas-phase species.
	void init(const Sprog::SpeciesPtrVector &sp);


	unsigned int m_i_AL;
	unsigned int m_i_AL2O3;
	unsigned int m_i_O2;
	
};
} // Processes
} // Sweep

#endif // SWP_ALUMINUM_SURFACE_OXIDATION_H_