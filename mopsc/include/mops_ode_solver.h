/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ODE_Solver class wraps the CVODE ODE solver to solve single
    reactor models.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#ifndef MOPS_ODE_SOLVER_H
#define MOPS_ODE_SOLVER_H

#include "mops_params.h"
#include "mops_src_terms.h"
#include "mops_reactor.h"
#include "mops_gpc_sensitivity.h"


// CVODE includes.
#include <sundials/sundials_dense.h>
#include "nvector/nvector_serial.h"
#include "cvodes_direct_impl.h"
#include "cvodes_impl.h" // For CVodeMem.
//#include "cvodes_dense_impl.h" // For DenseMat.

#include <istream>

#define ZERO  RCONST(0.0)

namespace Mops
{
class ODE_Solver
{
public:
    // Constructors.
    ODE_Solver(); // Default constructor.
    ODE_Solver(const ODE_Solver &copy);// Copy constructor.
    ODE_Solver(          // Stream-reading constructor.
        std::istream &in //   - Input stream.
        );

    // Destructor.
    ~ODE_Solver(void); // Default destructor.

    // Operators.
    // Assigned operator. This function cannot be used with sensitivity problem.
    // Sensitivity assignment has yet to be implemented.
    ODE_Solver &operator=(const ODE_Solver &rhs);

    // Enumeration of ODE solvers.
    // enum SolverType {CVODE_Solver, RADAU5_Solver};

    
    // SOLVER SETUP.

    // Initialises the solver at the given time.
    void Initialise(Reactor &reac);

    // Reset the solver.  Need to do this if the the reactor
    // contents has been changed between calls to Solve().
    void ResetSolver(void);

    // Reset the solver.  Need to do this if the the reactor
    // contents has been changed between calls to Solve().
    void ResetSolver(Reactor &reac);

    // Sets the time in the ODE solver.
//    void SetTime(double time);


    // RUNNING THE SOLVER.

    // Solves the reactor equations up to the given time, assuming
    // that it is in future to the current time.
	void Solve(Reactor &reac, double stop_time);


    // ERROR TOLERANCES.

    // Returns the absolute error tolerance used for ODE
    // calculations.
    double ATOL() const;

    // Sets the absolute error tolerance used for ODE
    // calculations.
    void SetATOL(double atol);

    // Returns the relative error tolerance used for ODE
    // calculations.
    double RTOL() const;

    // Sets the relative error tolerance used for ODE
    // calculations.
    void SetRTOL(double rtol);


    // EXTERNAL SOURCE TERMS.

    // Returns the vector of external source terms.
    const SrcProfile *const ExtSrcTerms(void) const;

    // Sets the external source terms.
    void SetExtSrcTerms(const SrcProfile &src);

    // Returns the external source term function.
    SrcTermFnPtr ExtSrcTermFn(void) const;

    // Sets the source term function pointer.
    void SetExtSrcTermFn(SrcTermFnPtr fn);


    // READ/WRITE/COPY FUNCTIONS.

    // Creates a copy of the solver object.
    ODE_Solver* Clone() const;

    // Writes the solver to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the solver data from a binary data stream.
    void Deserialize(std::istream &in);

    // RHS FUNCTION INTERFACE.
    // These function are implemented for being interface with Rhs functions of
    // CVODE and CVODES. Rhs functions are exported to mops_rhs_func.cpp/.h
    // so that the Rhs functions can be made non-static.

    // Return Reactor object.
    Reactor *GetReactor() { return m_reactor;}

    // Return a source term function pointer.
    SrcTermFnPtr GetsrcTermsFn() { return _srcTerms;};
    
    // Return a source term as SrcProfile object.
    const SrcProfile *GetsrcTerms() { return m_srcterms;};

    // Return number of equation.
    unsigned int GetNEquations() const { return m_neq; };

    //!Creates 2D array to store sensitivity values and initialises to zero.
    void InitialiseSensArray(int n_sensi, int n_species);

    //! Deletes the space allocated by InitialiseSensArray for the double** m_sensitivity.
    void DestroySensArray(int n_species);

    //!Transfers sensitivity values from m_yS into standard Sens 2D array.
    double** GetSensSolution(int n_sensi, int n_species);

    //! Returns number of sensitivities computed.
    unsigned int GetNSensitivities() const;

    // Set sensitivity object by making a copy of given sensitivity object.
    void SetSensitivity(Mops::SensitivityAnalyzer &sensi) const {m_sensi = sensi;};

    // Get sensitivity object.
    Mops::SensitivityAnalyzer &GetSensitivity() const {return m_sensi;};
	
private:
    //! Sensitivity matrix in double** for LOI access.
    double** m_sensitivity;

protected:
    // ODE solution variables.
    double m_rtol, m_atol;    // Relative and absolute tolerances.
    unsigned int m_neq;     // Number of equations solved.
//    unsigned int m_nsp;     // Number of species in current mechanism.
//    int m_iT;               // Index of temperature in solution vectors.
//    int m_iDens;            // Index of density in solution vectors.

    // Solution variables.
    double m_time;        // Current solution time.
    Reactor *m_reactor; // The reactor being solved.
    double *m_soln;       // Pointer to solution array (comes from Reactor object).
    double *m_deriv;      // Array to hold current solution derivatives.

    // External source terms.
    const SrcProfile *m_srcterms; // Vector of externally defined source terms on  the RHS.
    SrcTermFnPtr _srcTerms; // Source term function pointer.

    // Sensitivity related variables
    mutable Mops::SensitivityAnalyzer m_sensi;
    // Space required by rate parameter sensitivity.
    N_Vector *m_yS;


private:
    // CVODE variables.
    void *m_odewk;     // CVODE workspace.
    N_Vector m_solvec; // Internal solution array for CVODE interface.
    N_Vector m_yvec;   // Internal y work space for CVODE interface.


    // INITIALISATION AND DESTRUCTION.
    
    // Initialises the solver to the default state.
    void init(void);

    // Releases all memory used by the solver object.
    void releaseMemory(void);

    // Initialises the CVode ODE solver assuming that the
    // remainder of the the solver has been correctly set up.
    void InitCVode(void);
    
};
};

#endif
