//************************
// this file is to calculate the mass diffusion rate 
//and the influence of oxide cap.
#ifndef SWP_ALUMINUM_MASS_DIFFUSION_H_
#define SWP_ALUMINUM_MASS_DIFFUSION_H_

#include "swp_surface_reaction.h"

namespace Sweep {
namespace Processes {

class AlMassDiffusion: public Sweep::Processes::SurfaceReaction
{
public:
	// Mechanism constructor
	AlMassDiffusion(const Sweep::Mechanism &mech);

	// Copy constructor
	AlMassDiffusion(const AlMassDiffusion &copy);

	//stream reading constructor
	AlMassDiffusion(std::istream &in,
	const Sweep::Mechanism &mech);

	// Create a copy of the particle process.
	AlMassDiffusion *const AlMassDiffusion::Clone(void) const;

	// Assignment operator
	AlMassDiffusion &operator=(const AlMassDiffusion &rhs);

	// return the process type
	ProcessType ID(void) const;

	double CalcCoverFrac(const Cell &sys) const;

	// return the rate of the process for the given system.
	double Rate(
		double t,
		const Cell &sys,
		const Geometry::LocalGeometry1d &local_geom) const;
   
	// returns rate of the process for the given particle.
	double Rate(
		double t, 
		const Cell &sys,
		const Particle &sp) const;

	//writes the object to a binary stream.
	void Serialize(std::ostream &out) const;

	//read the object from a binary stream.
	void Deserialize(
		std::istream &in,
		const Sweep::Mechanism &mech
		);

private:
	// Private default constructor
	AlMassDiffusion();

	// Initialize indices for gas-phase species.
	void init(const Sprog::SpeciesPtrVector &sp);
	double m_fraction;
	// Index of Al and Al2O3 in the gas phase
	unsigned int m_i_AL;
	unsigned int m_i_AL2O3;

};
}
}
#endif // SWP_ALUMINUM_MASS_DIFFUSION_H_