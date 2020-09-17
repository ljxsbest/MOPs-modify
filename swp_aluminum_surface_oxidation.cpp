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
	m_i_O2(){}

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

// return the surface oxidation rate constant, m3.
double AlSurfOxidation::SurfaceRateConstant(double t,
										const Cell &sys) const
{  
	// Boltzmann constant, J/K.
	double kb = 1.38064825E-23; 
	// mass of single O2 molecule, kg. 
	double mass_O2 = 0.032 / 6.02214076E+23; 
	
	double V = sqrt(8 * kb * 
			sys.GasPhase().Temperature() / (mass_O2 * 3.141592653)); // m/s
	// O2 consume rate constant k1, m3.
	// where 0.1 is the ratio of reactant O2 and total O2. 
	// GetTotalDiameter2() = 4*r*r
	double k1= 0.1 * V * t 
				* PI * sys.Particles().GetTotalDiameter2(); 

	return k1;
}

// return O2 consume during t, mol
double AlSurfOxidation::O2ConsumeRate(double t,
										const Cell &sys) const
{
	double mole_O2Consume = SurfaceRateConstant(t, sys) * 
								sys.GasPhase().SpeciesConcentration(m_i_O2);

	return mole_O2Consume; // mol.
}

//  4/3 Al + O2 = 2/3 Al2O3.
double AlSurfOxidation::AlConsumeRate(double t,
										const Cell &sys) const
{
	double mole_AlConsume = 4 * SurfaceRateConstant(t, sys) * 
								sys.GasPhase().SpeciesConcentration(m_i_O2) / 3;
	return mole_AlConsume; // mol. 
}

double AlSurfOxidation::Al2O3ProduceRate(double t,
										const Cell &sys) const
{
	double mole_Al2O3Produce = mole_O2Consume * 4 / 3;
	return mole_Al2O3Produce; // mol.
}

//Now we get the O2 consume rate, (mol)
double AlSurfOxidation::Rate(double t,
							const Cell &sys) const
{
	double rate(1.0);
	// default reactor volume:
	double volume(1.0); // cm3.
	if (sys.GasPhase().Temperature() < 2700){
		// IUPAC's definition
		rate = O2ConsumeRate(t, sys) / (3 * volume * t);
	}
	else {
		rate = 0;
	}
	return rate;
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