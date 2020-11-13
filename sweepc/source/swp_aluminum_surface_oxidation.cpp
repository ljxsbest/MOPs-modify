#include "swp_aluminum_surface_oxidation.h"
#include "swp_mechanism.h"

namespace Sweep{
namespace Processes{

// all kinds of constructor
// Default constructor (private)
AlSurfOxidation::AlSurfOxidation()
:SurfaceReaction(),
	m_i_AL(0u),
	m_i_AL2O3(0u),
	m_i_O2(0u){}

//Mechanism Constructor
AlSurfOxidation::AlSurfOxidation(
	const Sweep::Mechanism &mech)
:SurfaceReaction(mech),
	m_i_AL(0u),
	m_i_AL2O3(0u),
	m_i_O2(0u)
{
	// call intialise to assign the indices
	init(*mech.Species());
}

// Copy constructor
AlSurfOxidation::AlSurfOxidation(const AlSurfOxidation &copy)
{
	*this = copy;
}

// Stream-reading constructor
AlSurfOxidation::AlSurfOxidation(
	std::istream &in,
	const Sweep::Mechanism &mech)
{
	Deserialize(in, mech);
}

// Assignment operator
AlSurfOxidation &AlSurfOxidation::operator =(const AlSurfOxidation &rhs)
{
	if (this != &rhs) {
		SurfaceReaction::operator =(rhs);
		m_i_AL = rhs.m_i_AL;
		m_i_AL2O3 = rhs.m_i_AL2O3;
		m_i_O2 = rhs.m_i_O2;
	}
	return *this;

}

// Create a copy of the particle process.
AlSurfOxidation *const AlSurfOxidation::Clone() const
{
	return new AlSurfOxidation(*this);
}

// return the process type
ProcessType AlSurfOxidation::ID() const { return AlSurfOxidation_ID; }

// assign species to object
void AlSurfOxidation::init(const Sprog::SpeciesPtrVector &sp)
{
	//Loop over species to find AL and AL2O3
	for (unsigned int i = 0; i !=sp.size(); ++i)
	{
		if ((sp[i]->Name() == "AL(l)") || (sp[i]->Name() == "AL(L)"))
			m_i_AL = i;
		if ((sp[i]->Name() == "AL2O3(l)") || (sp[i]->Name() == "AL2O3(L)")) 
			m_i_AL2O3 = i;
		if (sp[i]->Name() == "O2")
			m_i_O2 = i;
	}
	// Check if they are assigned
	if ((m_i_AL == 0u) && (m_i_AL2O3 == 0u) && (m_i_O2 == 0u))
		throw std::runtime_error("could not find AL, "
			"AL2O3 and O2 in gas phase in AlSurfOxidation::init");
}

// return the surface oxidation rate constant, m3/s.
double AlSurfOxidation::SurfaceRateConstant(double t,
										const Cell &sys,
										const Particle &sp) const
{  
	// Boltzmann constant, J/K. 
	// mass of single O2 molecule, kg. 
	double mass_O2 = 0.032 / NA;
	double v = sqrt(8 * KB * 
			sys.GasPhase().Temperature() / (mass_O2 * PI)); // m/s
	// O2 consume rate constant k1, m3.
	// where 0.1 is the ratio of reactant O2 and total O2. 
	// GetTotalDiameter2() = 4*r*r
	double k1 = 0.1 * v * sp.SurfaceArea(); 
	return k1;
}


// return O2 production/consumpation rate.
// default volume is 1.0, dimensionless.
double AlSurfOxidation::O2ConsumeRate(double t,
										const Cell &sys,
										const Particle &sp) const
{
	double mole_O2Consume = SurfaceRateConstant(t, sys, sp) *
		sys.GasPhase().SpeciesConcentration(m_i_O2) / sys.SampleVolume();
	return mole_O2Consume; // number/(m3*s).
}

//  4/3 Al + O2 = 2/3 Al2O3.
double AlSurfOxidation::AlConsumeRate(double t,
										const Cell &sys,
										const Particle &sp) const
{
	double mole_AlConsume = 4 * SurfaceRateConstant(t, sys, sp) * 
		sys.GasPhase().SpeciesConcentration(m_i_O2) / (3 * sp.Volume());
	return mole_AlConsume; //mol/s. 
}

double AlSurfOxidation::Al2O3ProduceRate(double t,
										const Cell &sys,
										const Particle &sp) const
{
	double mole_Al2O3Produce = 2 * SurfaceRateConstant(t, sys, sp) *
		sys.GasPhase().SpeciesConcentration(m_i_O2) / (3 * sp.Volume());
	return mole_Al2O3Produce; //mol/(m3*s).
}

//Now we get the O2 consume rate, (mol)
// based on former, rate = number/(m3*s)
double AlSurfOxidation::Rate(double t,
							const Cell &sys,
							const Particle &sp) const
{
	double rate(1.0);
	rate = O2ConsumeRate(t, sys, sp) / 3;
	// default reactor volume:
	if ((sys.GasPhase().Temperature() < 2710) && (sys.GasPhase().Temperature() >= 2350)){
		// IUPAC's definition
		rate *= 1;
	}
	else if (sys.GasPhase().Temperature() < 2350){
		rate *= 0.5;
	}
	else {
		rate = 0;
	}
	return rate;
}

// calculate the enthalpy:
// H/RT = a1 + a2*T/2 + a3 * pow(T, 2)/3 + a4 * pow(T, 3)/4 +a5 * pow(T, 4)/5 + a6/T 
// return molar production rate * enthalpy.
double AlSurfOxidation::HeatProdRate(double t,
					 const Cell &sys,
					 const Particle &sp) const
{
	double H_Al(1.0), H_O2(1.0), H_Al2O3(1.0);
	double rate_heatprod(1.0);
	H_Al = 3.83089866E+00 - 2.09027129E-05 * sys.GasPhase().Temperature()/2 + 
		1.04271684E-08 * pow(sys.GasPhase().Temperature(), 2) / 3 -
		2.04841051E-12 * pow(sys.GasPhase().Temperature(), 3) / 4 +
		1.39565517E-16 * pow(sys.GasPhase().Temperature(), 4) / 5 -
		9.97961566E+01 / sys.GasPhase().Temperature();

	H_Al *= R * sys.GasPhase().Temperature();

	H_O2 = 3.45852381E+00 + 1.04045351E-03 * sys.GasPhase().Temperature()/2 - 
		2.79664041E-07 * pow(sys.GasPhase().Temperature(), 2) / 3 +
		3.11439672E-11 * pow(sys.GasPhase().Temperature(), 3) / 4 -
		8.55656058E-16 * pow(sys.GasPhase().Temperature(), 4) / 5 +
		1.02229063E+04 / sys.GasPhase().Temperature();

	H_O2 *= R * sys.GasPhase().Temperature();

	H_Al2O3 = 1.95922550E+01 - 1.02229063E+04 / sys.GasPhase().Temperature();
	H_Al2O3 *= R * sys.GasPhase().Temperature();
    
    // calculate the heat production rate, J/s
    rate_heatprod = H_Al2O3 * Al2O3ProduceRate(t, sys, sp) -
    					H_Al * AlConsumeRate(t, sys, sp) -
    					H_O2 * O2ConsumeRate(t, sys, sp);

    return rate_heatprod;
}

double AlSurfOxidation::SetHeat(double a){
	return a;
}

//Write the object to a binary stream.
//@param out
void AlSurfOxidation::Serialize(std::ostream &out) const
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

		val = m_i_O2;
		out.write((char*)&val, sizeof(val));
	} else {
		throw invalid_argument("Output stream not ready"
			"in AlSurfOxidation::Serialize");
	}
}
//Read objects from a binary stream
//@param in and mech
void AlSurfOxidation::Deserialize(
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

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_i_O2 = (unsigned int) val;        

    } else {
        throw invalid_argument("Input stream not ready in AlSurfOxidation::Deserialize");
	}
}

}  //Processes
}  //Sweep