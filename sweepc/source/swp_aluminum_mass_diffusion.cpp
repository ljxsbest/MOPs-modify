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
		if (sp[i]->Name() == "AL2O3(L)") 
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
	//double r_cap(1.0);
	//double r_Al(1.0);
	//volume of Al and cap particle.
	double V_Al = sys.GasPhase().SpeciesConcentration(m_i_AL) /
		(2.377 * 27.0);
	double V_cap = sys.GasPhase().SpeciesConcentration(m_i_AL2O3) /
		(3.9 * 102.0);
	double r_Al = pow(3 * V_Al / (4 * PI), 1.0/3);
	double r_cap = pow(3 * V_cap / (4 * PI), 1.0/3);
	// the distance of centre of mass between two particles. 
	double m_fraction(1.0);
	if (r_Al >= r_cap){
		m_fraction = pow(r_cap / (2 * r_Al), 2); //should we define "georatio" here??
	}
	else if (r_Al < r_cap){
		m_fraction = 0.5 - r_Al / (4 * r_cap);
	}
	else{
		cout << "ErrGeo";
	}
	printf("%f\n ", m_fraction);
	return m_fraction;
}

// returns the rate of the process for the single particle. 
//（single particle?）
double AlMassDiffusion::Rate(double t, const Cell &sys, const Particle &sp) const
{
	double Alpha = 1;
	// surface energy from Storozhev2013.
	double SE = 0.7; // J/m2.
	// mass of Al molecule 
	double mass_Al = 0.027 / NA;
	// radius of Al particle, SA = 4 * pi * r * r.
	double r_Al = sqrt(sp.SurfaceArea() / (4 * PI));
	// velocity of Al molecule
	double v_Al = sqrt(8 * KB *
		sys.GasPhase().Temperature() / (mass_Al * PI));
	// density of Al particle.  kg / m3
	double Density = 2377;
	// equilibrium vapor pressure
	double p_e_r = sys.GasPhase().Pressure() *
		exp(2 * SE * mass_Al / (r_Al * Sweep::KB * 2377 * sys.GasPhase().Temperature()));
	// rate: number of Al atom / s
	double rate = Alpha * p_e_r * v_Al * sp.SurfaceArea() / (4 * KB * sys.GasPhase().Temperature());
	// rate: number/(m3 * s)
	// we can get different temeprature by choose 
	// Sp volume or sample volume.
	rate *= 1 / (NA * sys.SampleVolume());
	// geometric model
	rate *= 1 - CalcCoverFrac(sys);
	//choose mechanism based on temperature.
	if (sys.GasPhase().Temperature() >= 2690){
		rate *= 1;
	}
	else {
		rate = 0;
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