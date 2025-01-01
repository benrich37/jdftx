/*-------------------------------------------------------------------
Copyright 2014 Deniz Gunceler, Ravishankar Sundararaman

This file is part of JDFTx.

JDFTx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JDFTx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with JDFTx.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#ifndef JDFTX_ELECTRONIC_DUMP_INTERNAL_H
#define JDFTX_ELECTRONIC_DUMP_INTERNAL_H

#include <core/ScalarFieldArray.h>
#include <core/Coulomb.h>

class Everything;
class ColumnBundle;

//! @addtogroup Output
//! @{
//! @file Dump_internal.h Implementation internals for output modules


//-------------------- Implemented in DumpSIC.cpp ---------------------------

//! Output self-interaction correction for the KS eigenvalues
class DumpSelfInteractionCorrection
{
public:
	DumpSelfInteractionCorrection(const Everything& e);
	double operator()(std::vector<diagMatrix>* correctedEigenvalues); //!< Evaluates the self interaction energy and (optionally) returns the corrected band eigenvalues
	void dump(const char* filenamePattern);
	bool needsTau;  //!< The kinetic energy density is needed for meta-gga functionals.
private:
	const Everything& e;
	double calcSelfInteractionError(int q, int n); //!< Calculates the self-interaction error of the KS orbital atthe n'th band at q'th quantum number
	std::vector<ColumnBundle> DC; //!< ColumnBundle for the derivative of the wavefunctions in each cartesian direction
};

//---------------- Implemented in DumpExcitationsMoments.cpp -----------------

//! Dump information about excitation energies and matrix elements
void dumpExcitations(const Everything& e, const char* filename);

//! Dump coulomb matrix elements in FCIDUMP format
void dumpFCI(const Everything& e, const char* filename);

//! Dump dipole moments
void dumpMoment(const Everything& e, const char* filename);

namespace XC_Analysis
{	ScalarFieldArray tauWeizsacker(const Everything& e); //!< Output Weizsacker KE density
	ScalarFieldArray spness(const Everything& e); //!< Output 'single-particle-ness'
	ScalarFieldArray sHartree(const Everything& e); //!< output spin Hartree potentials
}

//! Dump band projections to atomic orbitals or ortho-orbitals depending on ortho, and complex or real based on norm
void dumpProjections(const Everything& e, const char* filename, bool ortho, bool norm);

//! Dump projections between atomic orbitals or ortho-orbitals depending on ortho, and complex or real based on norm
void dumpProjectionOverlap(const Everything& e, const char* filename, bool ortho, bool norm);

//---------------- Implemented in DumpChargedDefects.cpp -----------------

//! Slab dielectric function calculator
struct SlabEpsilon
{	string dtotFname; //!< reference electrostatic potential filename
	double sigma; //!< smoothing width
	vector3<> Efield; //!< reference electric field
	
	void dump(const Everything& e, ScalarField d_tot) const;
};

//! Bulk dielectric function calculator
struct BulkEpsilon
{	string dtotFname; //!< reference electrostatic potential filename
	vector3<> Efield; //!< reference electric field
	
	void dump(const Everything& e, ScalarField d_tot) const;
};

//! Charged defect correction calculator
struct ChargedDefect
{	//! Model charge decsription (one for each defect in unit cell)
	struct Center
	{	vector3<> pos; //!< defect center in lattice coordinates
		double q; //!< defect electron-count
		double sigma; //!< model charge Gaussian width
	};
	std::vector<Center> center; //!< list of defect positions in unit cell
	
	CoulombParams::Geometry geometry; //!< geometry of Coulomb interaction used for correction (could differ from main calculation)
	int iDir; //!< slab trunctaion direction
	
	string dtotFname; //!< electrostatic potential from reference neutral calculation
	
	double bulkEps; //!< bulk dielectric constant (Bulk mode only)
	string slabEpsFname; //!< slab dielectric profile (Slab mode only)
	
	double rMin; //!< Minimum distance from defect used for calculating alignment
	double rSigma; //!< Turn-on width of region used for calculating alignment
	
	void dump(const Everything& e, ScalarField d_tot) const;
};

//-------------------- Implemented in DumpCprime.cpp ---------------------------

struct DumpCprime
{
	double dk;
	double degeneracyThreshold;
	double vThreshold;
	bool realSpaceTruncated;
	
	DumpCprime(double dk=1E-4, double degeneracyThreshold=1E-6, double vThreshold=1E-4, bool realSpaceTruncated=true);
	void dump(Everything& e) const;

private:
	ColumnBundle getCprime(Everything& e, int q, int iDir, matrix& CprimeOC) const;
	ColumnBundle getCpert(Everything& e, int q, vector3<> dkVec, const matrix& dkDotV, matrix& CpertOC) const;
	matrix fixUnitary(const matrix& CpertOC, const diagMatrix& E, const matrix& dH) const;
};

//! @}
#endif // JDFTX_ELECTRONIC_DUMP_INTERNAL_H
