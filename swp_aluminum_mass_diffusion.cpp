#include "swp_aluminum_mass_diffusion.h"
#include "swp_mechanism.h"

namespace Sweep {
namespace Processes {

// Default constructor (private)
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
		if ((sp[i]->Name() == "AL(l)") || (sp[i]->Name() == "AL(L)"))
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

// return the process type
ProcessType AlMassDiffusion::ID() const { return AluminumSR_ID; }

// calculate the coverage fraction of Oxide cap.Casper
// we assume that the centre of cap is on the surface of Al particle.
double AlMassDiffusion::CalcCoverFrac(const Cell &sys) const
{
	//the molar mass of Al and Al2O3: 27 and 102.
	//the density of liquid Al and alpha Al2O3 are 2.377 g/cm3 and 3.9 g/cm3.
	double r_cap(1.0);
	double r_Al(1.0);
	//volume of Al and cap particle.
	double V_Al = sys.GasPhase().SpeciesConcentration(m_i_AL) /
		(2.377 * 27);
	double V_cap = sys.GasPhase().SpeciesConcentration(m_i_AL2O3) /
		(3.9 * 102);
	r_Al = pow(3 * V_Al / (4 * PI), 1 / 3);
	r_cap = pow(3 * V_cap / (4 * PI), 1 / 3);
	// the distance of centre of mass between two particles. 
	double gap = 1 - 0.5 * pow(r_cap / r_Al, 2);
	// free surface area of Al particle
	double Area_Al = 2 * PI * r_Al * (r_Al + gap);
	// 1 - the ratio of free SA and total SA. 
	double m_fraction = 1 - Area_Al / (4 * PI * pow(r_Al, 2));
	if (m_fraction > 0.5){
		m_fraction = 0.5; //should we define "georatio" here??
	}
	else if (m_fraction <= 0.5){
		m_fraction *= 1;
	}
	else{
		cout << "error Georetio";
	}
	return m_fraction;
}

// return the rate of the process for the given system.
double AlMassDiffusion::Rate(
	double t,
	const Cell &sys,
	const Geometry::LocalGeometry1d &local_geom) const
{
	double rate(1.0);

	//double T = gas.Temperature();
	double k1 = 1.0e8 ; //?????? needed to be verified.
	rate = k1*sys.GasPhase().SpeciesConcentration(m_i_AL);
	rate *= 1 - CalcCoverFrac(sys);
	//choose mechanism based on temperature
	if (sys.GasPhase().Temperature() >= 2700){
		rate *= 1;
	}
	else {
		rate = 0;
	}
	return rate;
}


// returns the rate of the process for the single particle. 
//（single particle?）
double AlMassDiffusion::Rate(double t, const Cell &sys, const Particle &sp) const
{
	double rate(1.0); 

	//double T = gas.Temperature();
	double k1 = 1.0e8; //?????? needed to be verified.
	rate = k1*sys.GasPhase().SpeciesConcentration(m_i_AL);
	rate *= 1 - CalcCoverFrac(sys);
	//choose mechanism based on temperature
	if (sys.GasPhase().Temperature() >= 2700){
		rate *= 1;
	}
	else {
		rate *= 0.1;
	}
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