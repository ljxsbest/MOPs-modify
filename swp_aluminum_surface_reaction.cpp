#include "swp_aluminum_surface_reaction.h"
#include "swp_mechanism.h"

namespace Sweep {
namespace Processes {

// Default constructor (private): 初始化列表？？？
AlMassDiffusion::AlMassDiffusion()
:SurfaceReaction(),
	m_i_AL(0u),
	m_i_AL2O3(0u){}


//Mechanism Constructor
AlMassDiffusion::AlMassDiffusion(
	const Sweep::Mechanism &mech)
:SurfaceReaction(mech),
	m_i_AL(0u),
	m_i_AL2O3(0u)
{
	// call intialise to assign the indices
	init(*mech.Species());
}

//Assigns gas-phase indices to object
void AlMassDiffusion::init(const Sprog::SpeciesPtrVector &sp)
{
	//Loop over species to find AL and AL2O3
	for (unsigned int i = 0; i !=sp.size(); ++i)
	{
		if (sp[i]->Name() == "AL")
			m_i_AL = i;
		if (sp[i]->Name() == "AL2O3")
			m_i_AL2O3 = i;
	}

	// Check if they are assigned
	if ((m_i_AL == 0u) && (m_i_AL2O3 == 0u))
		throw std::runtime_error("could not find AL and AL2O3 in gas phase in AlMassDiffusion::init");
}

// Copy constructor
AlMassDiffusion::AlMassDiffusion(const AlMassDiffusion &copy)
{
	*this = copy;
}

// Stream-reading constructor
AlMassDiffusion::AlMassDiffusion(
	std::istream &in,
	const Sweep::Mechanism &mech)
{
	Deserialize(in, mech);
}

// Assignment operator
AlMassDiffusion &AlMassDiffusion::operator =(const AlMassDiffusion &rhs)
{
	if (this != &rhs) {
		SurfaceReaction::operator =(rhs);
		m_i_AL = rhs.m_i_AL;
		m_i_AL2O3 = rhs.m_i_AL2O3;
	}
	return *this;

}

// Create a copy of the particle process.
AlMassDiffusion *const AlMassDiffusion::Clone() const
{
	return new AlMassDiffusion(*this);
}

//Return the process type
ProcessType AlMassDiffusion::ID() const {return AluminumSR_ID;}

//returns the rate of the process for the given system 
//（single particle?）
double AlMassDiffusion::Rate(double t, const Cell &sys, const Particle &sp) const
{
	double rate(1.0); //

	//double T = gas.Temperature();
	double k1 = 1.0e14; // constant diffusion rate.
	rate = k1*sqrt(sys.GasPhase().SpeciesConcentration(m_i_AL));
	return rate;
}

//Write the object to a binary stream.
//@param out
void AlMassDiffusion::Serialize(std::ostream &out) const
{
	if (out.good()) {
		//serialize base class.
		SurfaceReaction::Serialize(out);

		//write derived class data
		unsigned int val(0);

		val = m_i_AL;
		out.write((char*)&val, sizeof(val));

		val = m_i_AL2O3;
		out.write((char*)&val, sizeof(val));
	} else {
		throw invalid_argument("Output stream not ready"
			"in AlMassDiffusion::Serialize");
	}
}

//Read objects from a binary stream
//@param in and mech
void AlMassDiffusion::Deserialize(
	std::istream &in,
	const Sweep::Mechanism &mech)
{
	if (in.good()) {
		// Deserialize base class
		SurfaceReaction::Deserialize(in, mech);

		// Read derived class data
		unsigned int val(0);

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_i_AL = (unsigned int) val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_i_AL2O3 = (unsigned int) val;

    } else {
        throw invalid_argument("Input stream not ready in AlMassDiffusion::Deserialize");
	}
}

	
} //Processes

} //Sweep