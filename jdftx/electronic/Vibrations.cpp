/*-------------------------------------------------------------------
Copyright 2013 Ravishankar Sundararaman

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

#include <electronic/Vibrations.h>
#include <electronic/IonicMinimizer.h>
#include <electronic/Everything.h>
#include <core/LatticeUtils.h>
#include <core/Units.h>


Vibrations::Vibrations() : dr(0.01), centralDiff(false), useConstraints(false),
translationSym(true), rotationSym(false), omegaMin(2e-4), T(298*Kelvin), omegaResolution(1e-4), 
dumpK(false), collectConfigurations(true), iConfiguration(-1), nConfigurations(-1),
iConfigStart(0), iConfigStop(0), computeOnly(false)
{
}

void Vibrations::setup(Everything* e)
{	this->e = e;
	//Perform any compatibility checks here (so that dry runs will pick these up)
	if((translationSym || rotationSym) && useConstraints)
	{	logPrintf("WARNING: Vibrations: switching off translationSym and rotationSym since useConstraints is on.");
		translationSym = false;
		rotationSym = false;
	}
}

inline void setPtest(size_t iStart, size_t iStop, const vector3<int>& S, std::vector<double*> Ptest, vector3<> split)
{	vector3<> invS; for(int k=0; k<3; k++) invS[k] = 1./S[k];
	THREAD_rLoop
	(	for(int k=0; k<3; k++)
		{	double xk = invS[k] * iv[k];
			if(xk > split[k]) xk-=1.;
			Ptest[k][i] = xk;
		}
	)
}


void Vibrations::calculate()
{
	logPrintf("------ Vibrations::calculate() -------\n");
	logPrintf("WARNING: Vibrations module is experimental. Please report bugs!\n");
	logPrintf("Compare results with and without symmetries and report discrepancies.\n");

	//Initialize data structure to hold all relevant information for the calculation:
	VibrationsData data;

	//Construct the modes to evaluate
	construct_modes(data);

	// Determine which configurations to compute, and whether to only compute or also collect into Hessian/dipole derivative matrices:

	//                      ||      nConfigurations >= 0           ||  nConfigurations == -1  
	//                      ||                                     ||        (default)
	// -------------------- || --------------------------------    || -------------------------
	//                      ||            Run iConfigs             ||
	// iConfiguration >= 0  ||      iConfiguration through         || Run iConfiguration only     
	//                      ||  iConfiguration+nConfigurations-1   ||                                 
	// -------------------- || --------------------------------    || -------------------------
	// iConfiguration == -1 ||    Run first nConfigurations        || Run all configurations          
	//      (default)       ||                                     ||        (default)                       
        
	logPrintf("iConfiguration = %d, nConfigurations = %d\n", iConfiguration, nConfigurations);
	bool startSet = (iConfiguration>=0);
	bool lenSet = (nConfigurations>=0);
	iConfigStart = startSet ? iConfiguration : 0;
	iConfigStop = lenSet ? iConfigStart + nConfigurations : data.nModes + 1;
	logPrintf("iConfigStart = %d, iConfigStop = %d\n", iConfigStart, iConfigStop);
	computeOnly = startSet or lenSet; //If either of these are set, do not bother in anything besides computing and writing grad and Pel
	logPrintf("computeOnly = %d\n", computeOnly);
	
	//Determine number of degrees of freedom:
	logPrintf("Degrees of freedom: %d total, %d symmetry-independent.\n", data.nModes, data.nPrimary);
	if(!data.nModes)
	{	//Exit, but produce output in the same format as if there were modes:
		logPrintf("0 imaginary modes, 0 modes within cutoff, 0 real modes.\n");
		logPrintf("\nVibrational free energy components at T = %lg K:\n", T/Kelvin);
		logPrintf("\tZPE:   %15.6lf\n", 0.);
		logPrintf("\tEvib:  %15.6lf\n", 0.);
		logPrintf("\tTSvib: %15.6lf\n", 0.);
		logPrintf("\tAvib:  %15.6lf\n", 0.);
		logPrintf("\n");
		return;
	}
	
	//Find inverse of each symmetry matrix:
	set_iRotInv(data);
	
	//Initialize dipole measuring vector field
	nullToZero(Ptest, e->gInfo);
	threadLaunch(setPtest, e->gInfo.nr, e->gInfo.S, Ptest.data(), getSplit());

	//Construct maps for mapping iConfig to perturbation d (ds), and iMode to configuration indices (iModeToConfigs):
	logPrintf("Constructing maps ...\n"); logFlush();
	construct_maps(data);
	logPrintf("Constructing maps - done.\n"); logFlush();

	//Initialize Hessian, dipole derivative matrix, and multiplicity vector:
	logPrintf("Initializing K/dP/mult.\n"); logFlush();
	data.K = zeroes(data.nModes, data.nModes); //force matrix
	data.dP = zeroes(data.nModes, 3); //dipole derivative
	data.mult = diagMatrix(data.nModes, 0.); //multiplicity in entries due to symmetrization
	{	logPrintf("imin/grad0/Pel0.\n"); logFlush();
		IonicMinimizer imin(*e);
		IonicGradient grad0;
		vector3<> Pel0;
		//Evaluate forces and dipole moment at equilibrium configuration (iConfig=0)
		logPrintf("Computing K/dP/mult ...\n"); logFlush();
		compute_or_collect_iConfig(imin, data.ds, grad0, Pel0, 0);
		logPrintf("Computing K/dP/mult - done.\n"); logFlush();
		logPrintf("Looping through modes/configurations ...\n"); logFlush();
		for(int iMode=0; iMode<data.nModes; iMode++){
			//Evaluate forces and dipole moment at configuration(s) corresponding to this mode, and add contributions to K and dP
			process_mode(imin, iMode, data, grad0, Pel0);
		}
		logPrintf(" ... done.\n"); logFlush();
		if(computeOnly){
			logPrintf("All requested configurations computed and written to disk. Exiting.\n");
			return;
		}
	}
	logPrintf("Clearing mult/ds from memory.\n"); logFlush();
	data.mult.clear();
	data.ds.clear();

	//Correct for multiple countings:
	logPrintf("Accounting for multiplicity.\n"); logFlush();
	account_for_multiplicity(data);
	//Fill in modes set by translation symmetry, if any:
	logPrintf("Filling in modes set by translation symmetry.\n"); logFlush();
	fill_in_trans_sym_modes(data);

	//Symmetrize force matrix:
	logPrintf("Symmetrizing force matrix.\n"); logFlush();
	logPrintf("\nRelative symmetry error in force matrix = %lg\n", 0.5*nrm2(data.K - dagger(data.K))/nrm2(data.K));
	data.K = dagger_symmetrize(data.K);

	//Project out translation / rotation modes:
	logPrintf("Projecting out tr/rot modes.\n"); logFlush();
	apply_projections(data);
	//Construct mass-weighted frequency-squared matrix and solve for eigenvalues and eigenvectors:
	logPrintf("Solving mass-weighted eigenmodes.\n"); logFlush();
	solve_modes(data);
	//Enumerate imaginary, zero, and real modes:
	logPrintf("Classifying eigenmodes.\n"); logFlush();
	process_solved_modes(data);
	//Print mode information to out file:
	logPrintf("Printing out eigenmode information.\n"); logFlush();
	print_modes(data);
	//Print vibrational free energy components to out file:
	logPrintf("Printing out free energy.\n"); logFlush();
	print_free_energy(data);
	//Dump Hessian if requested:
	if(dumpK)
	{	string fname = e->dump.getFilename("K");
		logPrintf("\nWriting force matrix K to '%s' ... \n", fname.c_str()); logFlush();
		FILE* fp = fopen(fname.c_str(), "wb");
		if(!fp) die("Error opening file for writing.\n");
		data.K.write(fp);
		fclose(fp);
	}	
}

void Vibrations::construct_modes(VibrationsData& data)
{	//Create a non-constraint which simplifies the logic below
	SpeciesInfo::Constraint nullConstraint;
	nullConstraint.moveScale = 1.;
	nullConstraint.type = SpeciesInfo::Constraint::None;
	const auto& species = e->iInfo.species;
	// const auto& sym = *(data.sym);
	const std::vector<SpaceGroupOp>& sym = e->symmUnperturbed.getMatrices();
	// const auto& atomMap = *(data.atomMap);
	const auto& atomMap = e->symmUnperturbed.getAtomMap();
	data.nPrimary = 0;
	data.foundTranslatable = false;
	for(unsigned s=0; s<species.size(); s++)
	{	const SpeciesInfo& sp = *(species[s]);
		std::vector<bool> isPrimary(sp.atpos.size(), true);  //whether atom is the first of a set related by symmetries
		for(unsigned a=0; a<sp.atpos.size(); a++)
		{	//Find allowed directions of motion for this atom (accounting for move flags and constraints)
			std::vector< vector3<> > nSet; //list of allowed cartesian directions
			const SpeciesInfo::Constraint& c = useConstraints ? sp.constraints[a] : nullConstraint;
			if(!c.moveScale) continue;
			switch(c.type)
			{	case SpeciesInfo::Constraint::None:
				{	nSet.push_back(vector3<>(1,0,0));
					nSet.push_back(vector3<>(0,1,0));
					nSet.push_back(vector3<>(0,0,1));
					break;
				}
				case SpeciesInfo::Constraint::Planar:
				{	vector3<> a = fabs(c.d[0])<fabs(c.d[1]) ? vector3<>(1,0,0) : vector3<>(0,0,1); //=> a cannot be parallel to d
					a -= c.d * (dot(c.d,a)/c.d.length_squared()); //now a is perpendicular to d
					vector3<> b = cross(c.d, a); //b is perpendicular to a and d
					nSet.push_back(a / a.length());
					nSet.push_back(b / b.length());
					break;
				}
				case SpeciesInfo::Constraint::Linear:
				{	nSet.push_back(c.d / c.d.length());
					break;
				}
				case SpeciesInfo::Constraint::HyperPlane:
				{	die("Hyperplane constraint not yet supported in vibrations.\n")
				}
			}
			//Check whether to fillin with translation symmetries:
			bool fromTranslation = false;
			if(translationSym && !data.foundTranslatable)
			{	bool singleton = true; //invariant under symmetries
				for(unsigned iRot=0; iRot<sym.size(); iRot++)
					if(atomMap[s][a][iRot] != int(a))
					{	singleton = false;
						break;
					}
				if(singleton)
				{	fromTranslation = true;
					data.foundTranslatable = true;
					isPrimary[a] = false;
				}
			}
			//Add modes for each of these directions:
			for(const vector3<>& n: nSet)
			{	Mode mode = { s, a, n, isPrimary[a], fromTranslation };
				data.modes.push_back(mode);
			}
			if(isPrimary[a]) data.nPrimary += nSet.size();
			//Unset primary flags of symmetric images:
			for(unsigned iRot=0; iRot<sym.size(); iRot++)
				isPrimary[atomMap[s][a][iRot]] = false;
		}
	}
	// Cast to a signed integer for ease of use later on
	int nModes = data.modes.size();
	data.nModes = nModes;
}



void Vibrations::construct_maps(VibrationsData& data)
{	logPrintf("Constructing maps - creating references.\n"); logFlush();
	int nConfigs = 1 + data.nPrimary * (centralDiff ? 2 : 1);
	data.iModeToConfigs.resize(data.nModes);
	data.ds.resize(nConfigs);
	const std::vector<Mode>& modes = data.modes;
	std::vector<std::vector<int>>& iModeToConfigs = data.iModeToConfigs;
	std::vector<IonicGradient>& ds = data.ds;
	//Set perturbation for equilibrium configuration (iConfig=0) to zero:
	IonicGradient& d = ds[0];
	d.init(e->iInfo);
	//Collect the displacements corresponding to each configuration
	int iConfig = 1;
	// int nModes = modes.size();
	logPrintf("Constructing maps - beginning loop.\n"); logFlush();
	for(int iMode=0; iMode<data.nModes; iMode++)
	{	const Mode& mode = modes[iMode];
		std::vector<int>& iConfigs = iModeToConfigs[iMode];
		IonicGradient& d = ds[iConfig];
		d.init(e->iInfo);
		if(!mode.isPrimary) continue;
		logPrintf("Constructing maps - pushing back iConfig to iConfigs.\n"); logFlush();
		iConfigs.push_back(iConfig);
		logPrintf("Constructing maps - creating d.\n"); logFlush();
		logPrintf("Constructing maps - modifying d.\n"); logFlush();
		d[mode.s][mode.a] = mode.n;
		ds[iConfig] = d;
		if(centralDiff)
		{	iConfig++;
			IonicGradient& d2 = ds[iConfig];
			d2.init(e->iInfo);
			d2[mode.s][mode.a] = -mode.n;
			iConfigs.push_back(iConfig);
			ds[iConfig] = d2;
		}
		iConfig++;
	}
}

void Vibrations::set_iRotInv(VibrationsData& data)
{	//Find inverse of each symmetry matrix:
	const std::vector<SpaceGroupOp>& sym = e->symmUnperturbed.getMatrices();
	data.iRotInv = std::vector<unsigned>(sym.size());
	// std::vector<unsigned>& iRotInv = data.iRotInv;
	// std::vector<unsigned> iRotInv(sym.size());
	for(unsigned iRot1=0; iRot1<sym.size(); iRot1++)
		for(unsigned iRot2=iRot1; iRot2<sym.size(); iRot2++)
			if(sym[iRot1].rot * sym[iRot2].rot == matrix3<int>(1,1,1))
			{	data.iRotInv[iRot1] = iRot2;
				data.iRotInv[iRot2] = iRot1;
				continue;
			}
}

void Vibrations::compute_iConfig(IonicMinimizer& imin, IonicGradient& d, IonicGradient& grad, vector3<>& Pel)
{	//Compute forces and dipole derivatives for this configuration
	logPrintf("compute_iConfig - stepping.\n"); logFlush();
	imin.step(d, dr);
	logPrintf("compute_iConfig - computing grad.\n"); logFlush();
	imin.compute(&grad, 0);
	logPrintf("compute_iConfig - computing Pel.\n"); logFlush();
	Pel = getPel();
	logPrintf("compute_iConfig - unstepping.\n"); logFlush();
	imin.step(d, -dr); //restore to unperturbed configuration
}

void Vibrations::compute_or_collect_iConfig(IonicMinimizer& imin, std::vector<IonicGradient>& ds, IonicGradient& grad, vector3<>& Pel, int iConfig)
{	//Compute forces and dipole derivatives for this configuration
	logPrintf("compute_or_collect_iConfig - setting fnames.\n"); logFlush();
	string grad_fname = e->dump.getFilename(string(("grad" + std::to_string(iConfig)).c_str()));
	string Pel_fname = e->dump.getFilename(string(("Pel" + std::to_string(iConfig)).c_str()));
	bool perform_compute = false;
	std::vector<string> req_fnames = {grad_fname, Pel_fname};
	logPrintf("compute_or_collect_iConfig - check for file pre-existence.\n"); logFlush();
	for(const string& fname: req_fnames)
	{	FILE* fp = fopen(fname.c_str(), "rb");
		if(!fp){
			perform_compute = true;
			break;
		}
		else fclose(fp);
	}
	if(perform_compute){
		logPrintf("compute_or_collect_iConfig - setting up compute for iConfig.\n"); logFlush();
		// IonicGradient& d = ds[iConfig];
		IonicGradient& d = ds.at(iConfig);
		compute_iConfig(imin, d, grad, Pel);
		const int& nConfigs = ds.size();
		logPrintf("Completed %d of %d configurations.\n", iConfig+1, nConfigs);
		logPrintf("compute_or_collect_iConfig - writing grad.\n"); logFlush();
		grad.write(grad_fname.c_str());
		matrix Pel_mat = zeroes(3,1);
		for (int k=0; k<3; k++){
			Pel_mat(k,0) = Pel[k];
		}
		logPrintf("compute_or_collect_iConfig - writing Pel.\n"); logFlush();
		Pel_mat.write(Pel_fname.c_str());
		if(iConfig == 0){
			e->dump(DumpFreq_Ionic, iConfig);
		}
	}
	else if(!computeOnly){
		logPrintf("compute_or_collect_iConfig - reading grad.\n"); logFlush();
		grad.read(grad_fname.c_str());
		matrix Pel_mat = zeroes(3,1);
		// matrix3<> Pel_mat(1, 3, 3);
		logPrintf("compute_or_collect_iConfig - reading Pel.\n"); logFlush();
		Pel_mat.read(Pel_fname.c_str());
		complex *Pel_data = Pel_mat.data();
		for (int k=0; k<3; k++){
			Pel[k] += Pel_data[Pel_mat.index(k,0)].real();
		}
	}
}

void Vibrations::process_mode(IonicMinimizer& imin, int iMode, VibrationsData& data, IonicGradient& grad0, vector3<>& Pel0){
	logPrintf("Processing mode - selecting mode.\n"); logFlush();
	// Mode mode = data.modes[iMode];
	Mode mode = data.modes.at(iMode);
	if(mode.isPrimary){
		logPrintf("Processing mode - setting/initializing references.\n"); logFlush();
		IonicGradient gradPlus, gradMinus;
		vector3<> PelPlus, PelMinus;
		std::vector<IonicGradient>& ds = data.ds;
		std::vector<int> iConfigs = data.iModeToConfigs[iMode];
		IonicGradient Kcur; 
		Kcur.init(e->iInfo);
		vector3<> dPcur;
		// If this is a compute-only run, just compute every configuration within the range as gradPlus/PelPlus and move on
		if(computeOnly){
			logPrintf("Processing mode - compute only - looping through configs for mode.\n"); logFlush();
			for(int iConfig: iConfigs){
				if(iConfig >= iConfigStart && iConfig < iConfigStop){
					logPrintf("Processing mode - compute only - compute or collecting.\n"); logFlush();
					compute_or_collect_iConfig(imin, ds, gradPlus, PelPlus, iConfig);
				}
			}
		}
		else{
			logPrintf("Processing mode - compute and analyze - looping through configs for mode.\n"); logFlush();
			int iiConfig = 0;
			// compute_or_collect_iConfig(imin, ds, gradPlus, PelPlus, iConfigs[iiConfig++]);
			compute_or_collect_iConfig(imin, ds, gradPlus, PelPlus, iConfigs.at(iiConfig++));
			if(centralDiff){
				// compute_or_collect_iConfig(imin, ds, gradMinus, PelMinus, iConfigs[iiConfig++]);
				compute_or_collect_iConfig(imin, ds, gradMinus, PelMinus, iConfigs.at(iiConfig++));
				Kcur = (gradPlus - gradMinus) * (0.5/dr);
				dPcur = (PelPlus - PelMinus) * (0.5/dr);
			}
			else{
				Kcur = (gradPlus - grad0) * (1./dr);
				dPcur = (PelPlus - Pel0) * (1./dr);
			}
			logPrintf("Processing mode - compute and analyze - correcting for ionic contribution to dipole.\n"); logFlush();
			const auto& species = e->iInfo.species;
			// dPcur -= species[mode.s]->Z * mode.n; //ionic contribution to dipole derivative
			dPcur -= species.at(mode.s)->Z * mode.n; //ionic contribution to dipole derivative
			logPrintf("Processing mode - compute and analyze - collecting current contributions to K and dP.\n"); logFlush();
			collect_cur_contributions(data, Kcur, dPcur, iMode);
		}
	}
}


void Vibrations::collect_cur_contributions(VibrationsData& data, const IonicGradient& Kcur, const vector3<>& dPcur, int iMode){
	logPrintf("collect_cur_contributions - setting up references.\n"); logFlush();
	const std::vector<SpaceGroupOp>& sym = e->symmUnperturbed.getMatrices();
	const auto& atomMap = e->symmUnperturbed.getAtomMap();
	const auto& modes = data.modes;
	// Mode mode = modes[iMode];
	Mode mode = modes.at(iMode);
	const auto& iRotInv = data.iRotInv;
	std::vector<double>& mult = data.mult;
	matrix& K = data.K;
	matrix& dP = data.dP;
	logPrintf("collect_cur_contributions - looping through iRot.\n"); logFlush();
	for(unsigned iRot=0; iRot<sym.size(); iRot++){	
		matrix3<> rot = e->gInfo.R * sym[iRot].rot * inv(e->gInfo.R); //cartesian rotation matrix corresponding to symmetry
		// matrix3<> rot = e->gInfo.R * sym.at(iRot).rot * inv(e->gInfo.R); //cartesian rotation matrix corresponding to symmetry
		//Modes corresponding to displacement (first index of matrix):
		// unsigned a1 = atomMap[mode.s][mode.a][iRot];
		unsigned a1 = atomMap.at(mode.s).at(mode.a).at(iRot);
		vector3<> n1 = rot * mode.n;
		std::map<int,double> dModes;
		logPrintf("collect_cur_contributions - looping through modes (1).\n"); logFlush();
		for(int i1=0; i1<data.nModes; i1++){
			// if(modes[i1].s==mode.s && modes[i1].a==a1)
			const Mode& mode1 = modes.at(i1);
			if(mode1.s==mode.s && mode1.a==a1){	
				// double w = dot(n1, modes[i1].n); //projection weight
				double w = dot(n1, mode1.n); //projection weight
				if(fabs(w) < symmThreshold) continue;
				mult[i1] += w*w; //symmetry multiplicity
				//Loop over modes corresponding to force (second index of matrix):
				logPrintf("collect_cur_contributions - looping through modes (2).\n"); logFlush();
				for(int i2=0; i2<data.nModes; i2++)
				{	const Mode& mode2 = modes[i2];
					logPrintf("collect_cur_contributions - getting K i1 i2 contrib.\n"); logFlush();
					logPrintf("collect_cur_contributions - inverseRotation.\n"); logFlush();
					unsigned inverseRotation = iRotInv.at(iRot);
					logPrintf("collect_cur_contributions - a2.\n"); logFlush();
					unsigned a2 = atomMap.at(mode2.s).at(mode2.a).at(inverseRotation);
					// unsigned a2 = atomMap[mode2.s][mode2.a][iRotInv[iRot]]; //index of atom which upon rotation rot maps onto atom mode2.a
					
					logPrintf("collect_cur_contributions - i1 = %d, i2 = %d, a2 = %d, mode2.s = %d, mode2.a = %d, iRotInv[iRot] = %d\n", i1, i2, a2, mode2.s, mode2.a, iRotInv[iRot]); logFlush();
					// logPrintf("Kcur species=%d/%d, atom=%d/%d, rotation=%d/%d\n",
					// 			mode2.s, Kcur.size(),
					// 			a2, Kcur.at(mode2.s).size(),
					// 			inverseRotation, iRotInv.size());
					logPrintf("Kcur species=%d/", mode2.s); logFlush();
					logPrintf("%d, ", Kcur.size()); logFlush();
					logPrintf("atom=%d/", a2); logFlush();
					logPrintf("%d, ", Kcur.at(mode2.s).size()); logFlush();
					logPrintf("rotation=%d/", inverseRotation); logFlush();
					logPrintf("%d\n", iRotInv.size()); logFlush();
					logPrintf("collect_cur_contributions - getting K i1 i2 contrib - rotKcur\n"); logFlush();
					vector3<> rotKcur = rot * Kcur.at(mode2.s).at(a2);
					// vector3<> rotKcur = rot * Kcur[mode2.s][a2]; //rotated force derivative
					logPrintf("collect_cur_contributions - getting K i1 i2 contrib - part2\n"); logFlush();
					double contrib = w * dot(mode2.n, rotKcur);
					logPrintf("collect_cur_contributions - adding to K.\n"); logFlush();
					K.data()[K.index(i1,i2)] += contrib;
				}
				//Dipole derivatives:
				vector3<> rot_dPcur = rot * dPcur; //rotated dipole derivative
				for(int k=0; k<3; k++){
					logPrintf("collect_cur_contributions - adding to dP.\n"); logFlush();
					dP.data()[dP.index(i1,k)] += w * rot_dPcur[k];
				}
			}
		}
	}
}

void Vibrations::account_for_multiplicity(VibrationsData& data){
	logPrintf("account_for_multiplicity - inverting mult.\n"); logFlush();
	//Invert multiplicity matrixZero out  modes to be set by translational symmetry:
	for(int i=0; i<data.nModes; i++)
		// data.mult[i] = data.modes[i].fromTranslation ? 0. : 1./data.mult[i];
		data.mult.at(i) = data.modes.at(i).fromTranslation ? 0. : 1./data.mult.at(i);
	
	//Correct for multiple counting:
	logPrintf("account_for_multiplicity - correcting K and dP.\n"); logFlush();
	data.K = data.mult * data.K;
	data.dP = data.mult * data.dP;
}

void Vibrations::fill_in_trans_sym_modes(VibrationsData& data){
	logPrintf("fill_in_trans_sym_modes - setting references.\n"); logFlush();
	const std::vector<Mode>& modes = data.modes;
	matrix& K = data.K;
	matrix& dP = data.dP;
	logPrintf("fill_in_trans_sym_modes - looping through modes.\n"); logFlush();
	// for(int i1=0; i1<nModes; i1++) if(modes[i1].fromTranslation)
	for(int i1=0; i1<data.nModes; i1++) if(modes.at(i1).fromTranslation)
	{	//Create a uniform unit displacement of all atoms which moves current atom according to mode:
		matrix x(1, data.nModes); complex* xData = x.data();
		for(int i2=0; i2<data.nModes; i2++)
			// xData[x.index(0,i2)] = dot(modes[i1].n, modes[i2].n);
			xData[x.index(0,i2)] = dot(modes.at(i1).n, modes.at(i2).n);
		//A uniform displacement of all atoms should yield no net force
		//Except three rows of K are zero; set them so that the above becomes true.
		logPrintf("fill_in_trans_sym_modes - setting fill-ins for K/dP.\n"); logFlush();
		K.set(i1,i1+1, 0,data.nModes, -(x * K));
		//Simiarly polarization due to uniform displacement should be zero:
		dP.set(i1,i1+1, 0,3, -(x * dP));
	}
}

void Vibrations::apply_projections(VibrationsData& data){
	logPrintf("apply_projections - setting references.\n"); logFlush();
	const std::vector<Mode>& modes = data.modes;
	matrix& K = data.K;
	//Project out translation / rotation modes:
	matrix projector(data.nModes, 6); int nProjectors=0;
	complex* projData = projector.data();
	if(translationSym)
	{	for(int k=0; k<3; k++)
		{	
			vector3<> e(0,0,0); e[k]=1; //unit vector
			for(int i=0; i<data.nModes; i++)
				// projData[projector.index(i,nProjectors)] = dot(modes[i].n, e);
				projData[projector.index(i,nProjectors)] = dot(modes.at(i).n, e);
			nProjectors++;
		}
	}
	if(rotationSym)
	{	IonicGradient r = getCMcoords();
		//Compute inertia tensor:	
		matrix3<> I;
		for(unsigned s=0; s<e->iInfo.species.size(); s++)
		{	const SpeciesInfo& sp = *(e->iInfo.species[s]);
			for(const vector3<>& rAtom: r[s])
				I += sp.mass * (rAtom.length_squared()*matrix3<>(1,1,1) - outer(rAtom,rAtom));
		}
		//Get principal axis and moments:
		matrix Imat(3,3); for(int j=0; j<3; j++) for(int k=0; k<3; k++) Imat.set(j,k, I(j,k));
		matrix Ievecs; diagMatrix Ieigs; Imat.diagonalize(Ievecs, Ieigs);
		complex* IevecsData = Ievecs.data();
		for(int j=0; j<3; j++)
			if(Ieigs[j] > symmThreshold)
			{	double meanPhase, sigmaPhase, rmsImagErr;
				removePhase(3, IevecsData+Ievecs.index(0,j), meanPhase, sigmaPhase, rmsImagErr);
				vector3<> axis; for(int k=0; k<3; k++) axis[k] = IevecsData[Ievecs.index(k,j)].real();
				//Add rotational projectors for each axis with non-zero moment:
				for(int i=0; i<data.nModes; i++){
					const Mode& modei = modes[i];
					projData[projector.index(i,nProjectors)] = box(modei.n, axis, r[modei.s][modei.a]);
				}
				nProjectors++;
			}
	}
	if(nProjectors)
	{	projector = projector(0,data.nModes, 0,nProjectors); //discard empty columns
		projector = projector * invsqrt(dagger(projector)*projector); //orthonormalize
		matrix ppDag = projector * dagger(projector);
		matrix IminPpdag = eye(data.nModes) - ppDag;
		K = IminPpdag * K * IminPpdag;
		//dP -= ppDag * dP;
		logPrintf("Projected out %d rotation+translation modes\n", nProjectors);
	}
}

void Vibrations::solve_modes(VibrationsData& data){
	logPrintf("solve_modes - setting references.\n"); logFlush();
	const auto& species = e->iInfo.species;

	//Initialize mass matrix:
	logPrintf("solve_modes - Initialize mass matrix.\n"); logFlush();
	data.invsqrtM = diagMatrix(data.nModes);
	for(int i=0; i<data.nModes; i++)
		data.invsqrtM.at(i) = 1./sqrt(species[data.modes.at(i).s]->mass * amu);
	// for(int i=0; i<nModes; i++)
	// 	data.invsqrtM[i] = 1./sqrt(species[data.modes[i].s]->mass * amu);
	
	//Construct and diagonalize frequency-squared matrix:
	// data.omegaSq = zeroes(nModes, nModes);
	logPrintf("solve_modes - Constructing omegaAq.\n"); logFlush();
	data.omegaSq = data.invsqrtM * data.K * data.invsqrtM;
	data.omegaSqEigs = diagMatrix(data.nModes);
	data.omegaSqEvecs = zeroes(data.nModes, data.nModes);
	logPrintf("solve_modes - diagonalize omegaSq.\n"); logFlush();
	data.omegaSq.diagonalize(data.omegaSqEvecs, data.omegaSqEigs);
}

void Vibrations::process_solved_modes(VibrationsData& data){
	logPrintf("process_solved_modes - setting references.\n"); logFlush();

	//Determine number of modes of each type:
	data.iZeroStart=0, data.iRealStart=data.nModes;
	double omegaMinSq = omegaMin*omegaMin;
	double omegaSqEigPrev = -DBL_MAX;
	logPrintf("process_solved_modes - counting modes by type.\n"); logFlush();
	for(int i=0; i<data.nModes; i++)
	{	
		// double omegaSqEig = data.omegaSqEigs[i];
		double omegaSqEig = data.omegaSqEigs.at(i);
		if(omegaSqEigPrev<-omegaMinSq && omegaSqEig>-omegaMinSq) data.iZeroStart=i;
		if(omegaSqEigPrev<+omegaMinSq && omegaSqEig>+omegaMinSq) data.iRealStart=i;
		omegaSqEigPrev = omegaSqEig;
	}
	logPrintf("%d imaginary modes, %d modes within cutoff, %d real modes.\n", data.iZeroStart, data.iRealStart-data.iZeroStart, data.nModes-data.iRealStart);
	
	//Detect degeneracies:
	logPrintf("process_solved_modes - detecting degeneracies.\n"); logFlush();
	for(int i=1; i<data.nModes; i++)
		// if(fabs(sqrt(fabs(data.omegaSqEigs[i])) - sqrt(fabs(data.omegaSqEigs[i-1]))) > omegaResolution)
		if(fabs(sqrt(fabs(data.omegaSqEigs.at(i))) - sqrt(fabs(data.omegaSqEigs.at(i-1)))) > omegaResolution)
			data.iFreqChange.insert(i);
	data.iFreqChange.insert(0);
	data.iFreqChange.insert(data.iZeroStart);
	data.iFreqChange.insert(data.iRealStart);
	data.iFreqChange.insert(data.nModes);
}

void Vibrations::print_modes(VibrationsData& data){
	logPrintf("print_modes - setting references.\n"); logFlush();
	const double fineStructConst = 7.29735257e-3;
	std::vector<Mode>& modes = data.modes;
	diagMatrix& invsqrtM = data.invsqrtM;
	diagMatrix& omegaSqEigs = data.omegaSqEigs;
	matrix& omegaSqEvecs = data.omegaSqEvecs;
	int& iZeroStart = data.iZeroStart;
	int& iRealStart = data.iRealStart;
	std::set<int>& iFreqChange = data.iFreqChange;

	//Print modes:
	matrix dEvecs = invsqrtM * omegaSqEvecs; //displacements of the eigenvectors
	matrix Pevecs = dagger(dEvecs) * data.dP; // dipole moments of the eigenvectors
	diagMatrix PsqEvecs = diag(Pevecs * dagger(Pevecs)); //dipole intensity of the eigenvectors
	complex* dEvecsData = dEvecs.data();
	for(int i=0; i<data.nModes; i++)
	{	//Classify mode:
		string modeType; int iMode=0;
		if(i<iZeroStart) { iMode=i; modeType = "Imaginary"; }
		else if(i<iRealStart) { iMode=i-iZeroStart; modeType = "Zero"; }
		else { iMode=i-iRealStart; modeType = "Real"; }
		//Header:
		logPrintf("\n%s mode %d:\n", modeType.c_str(), iMode+1);
		//Frequency:
		double omega = sqrt(fabs(omegaSqEigs[i]));
		const char* omegaSuffix = omegaSqEigs[i]<0 ? "i" : "";
		logPrintf("Frequency: %.6lf%s Eh [ %.0lf%s cm^-1 ]\n", omega, omegaSuffix, omega/invcm, omegaSuffix);
		//Degeneracy:
		auto iterStop = std::upper_bound(iFreqChange.begin(), iFreqChange.end(), i);
		auto iterStart = iterStop; iterStart--;
		int degeneracyCount = (*iterStop) - (*iterStart);
		int degeneracyIndex = i - (*iterStart);
		logPrintf("Degeneracy: %d of %d\n", degeneracyIndex+1, degeneracyCount);
		//IR intensity:
		logPrintf("IR intensity: %.4f e^2/amu [ %.1f km/mol ]\n", PsqEvecs[i]*amu,
			PsqEvecs[i] * (M_PI/3.)*pow(fineStructConst,2) / (1e3*meter/mol));
		//Displacements:
		double meanPhase, sigmaPhase, rmsImagErr;
		removePhase(3, dEvecsData+dEvecs.index(0,i), meanPhase, sigmaPhase, rmsImagErr);
		IonicGradient d; d.init(e->iInfo);
		for(int j=0; j<data.nModes; j++)
			d[modes[j].s][modes[j].a] += modes[j].n * dEvecsData[dEvecs.index(j,i)].real();
		logPrintf("Displacements:\n");
		for(unsigned s=0; s<d.size(); s++)
		{	const SpeciesInfo& sp = *(e->iInfo.species[s]);
			for(const vector3<>& r: d[s])
				logPrintf("disp %s %19.15lf %19.15lf %19.15lf\n", sp.name.c_str(), r[0], r[1], r[2]);
		}
	}
}

void Vibrations::print_free_energy(VibrationsData& data){
	logPrintf("print_free_energy - setting references.\n"); logFlush();
	//Global quantities:
	double ZPE = 0., Evib = 0., Avib = 0.; //zero-point energy, average energy and free energy
	for(int i=data.iRealStart; i<data.nModes; i++)
	{	double omega = sqrt(data.omegaSqEigs[i]);
		double expMomegaByT = exp(-omega/T);
		ZPE += 0.5*omega;
		Evib += 0.5*omega + omega * expMomegaByT / (1.-expMomegaByT);
		Avib += 0.5*omega + T * log(1.-expMomegaByT);
	}
	double TSvib = Evib - Avib;
	logPrintf("\nVibrational free energy components at T = %lg K:\n", T/Kelvin);
	logPrintf("\tZPE:   %15.6lf\n", ZPE);
	logPrintf("\tEvib:  %15.6lf\n", Evib);
	logPrintf("\tTSvib: %15.6lf\n", TSvib);
	logPrintf("\tAvib:  %15.6lf\n", Avib);
}

vector3<> Vibrations::getSplit() const
{	//Collect lattice coordinates of all atoms in [0,1)
	std::vector<double> x[3];
	for(const auto& sp: e->iInfo.species)
		for(const vector3<>& at: sp->atpos)
			for(int k=0; k<3; k++)
				x[k].push_back(at[k] - floor(at[k]));
	//Find best periodicity split point along each lattice direction:
	vector3<> split;
	for(int k=0; k<3; k++)
	{	std::sort(x[k].begin(), x[k].end());
		//Determine interval lengths:
		std::vector<double> dx(x[k].size());
		for(unsigned i=0; i<x[k].size(); i++)
			dx[i] = x[k][i] - (i ? x[k][i-1] : x[k].back()-1.);
		//Set split point to the midpoint of longest interval:
		int iSplit = std::max_element(dx.begin(), dx.end()) - dx.begin();
		split[k] = x[k][iSplit] - 0.5*dx[iSplit];
		split[k] -= floor(split[k]); //map to [0,1)
	}
	return split;
}

IonicGradient Vibrations::getCMcoords() const
{	vector3<> split = getSplit();
	//Collect cartesian atom coordinates with wrapping consistent with the above determined split:
	IonicGradient r; r.init(e->iInfo);
	vector3<> rMsum; double Msum = 0.; //sums for determining CM
	for(unsigned s=0; s<e->iInfo.species.size(); s++)
	{	const SpeciesInfo& sp = *(e->iInfo.species[s]);
		for(unsigned a=0; a<sp.atpos.size(); a++)
		{	vector3<> xWrapped = sp.atpos[a];
			for(int k=0; k<3; k++)
			{	xWrapped[k] -= floor(xWrapped[k]); //map to [0,1)
				if(xWrapped[k] > split[k]) xWrapped[k] -= 1.; //split periodicity along determined point
			}
			r[s][a] = e->gInfo.R * xWrapped; //to cartesian coords
			rMsum += r[s][a] * sp.mass;
			Msum += sp.mass;
		}
	}
	vector3<> rCM = rMsum / Msum; //center of mass
	//Shift coordinates relative to center of mass:
	for(std::vector< vector3<> >& rSpecies: r)
		for(vector3<>& rAtom: rSpecies)
			rAtom -= rCM;
	return r;
}

vector3<> Vibrations::getPel() const
{	vector3<> Pel;
	for(int k=0; k<3; k++)
		Pel[k] = e->gInfo.dV * dot(Ptest[k], e->eVars.get_nTot());
	return e->gInfo.R * Pel; //convert to Cartesian coordinates
}
