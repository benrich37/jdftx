#include <electronic/IonicGaussianPotential.h>

double IonicGaussianPotential::energyAndGrad(const GridInfo& gInfo, std::vector<Atom>& atoms, matrix3<>* E_RRTptr) const
{
	double E = 0.; matrix3<> E_RRT;
	double mhalfInvSigmaSq = -0.5/(sigma*sigma);
	for(Atom& atom: atoms)
		if(atom.sp == iSpecies)
		{
			//Compute wrapped Cartesian position:
			vector3<> posWrapped = atom.pos - center; //shift to origin of Gaussian
			for(int iDir=0; iDir<3; iDir++)
				posWrapped[iDir] -= floor(0.5 + posWrapped[iDir]); //to [-0.5, 0.5)
			vector3<> r = gInfo.R * posWrapped;
			
			//Extract appropriate reduced coordinate:
			double rSq = 0.; vector3<> rSq_r;
			switch(geometry)
			{	case Spherical: { rSq = r.length_squared(); rSq_r = 2 * r; break; }
				case Cylindrical: { rSq = r[0]*r[0] + r[1]*r[1]; rSq_r = vector3<>(2*r[0], 2*r[1], 0.); break; }
				case Planar: { rSq = r[2]*r[2]; rSq_r = vector3<>(0., 0., 2*r[2]); break; }
			}
			
			//Accumulate energy, force and stress (if needed):
			double Ecur = U0 * exp(mhalfInvSigmaSq * rSq);
			vector3<> force = (-Ecur * mhalfInvSigmaSq) * rSq_r; //Cartesian force
			E += Ecur;
			atom.force += gInfo.RT * force; //convert to contravariant force
			if(E_RRTptr) E_RRT -= outer(r, force); //lattice gradient is negative of virial
		}
	if(E_RRTptr) *E_RRTptr += 0.5*(E_RRT + (~E_RRT)); //store symmetrized stress contribution
	return E;
}


double DynamicPlanarGaussianPotential::cum_force_z(const GridInfo& gInfo, std::vector<Atom>& atoms, matrix3<>* E_RRTptr, vector3<>* center) const
{
	double E = 0.; matrix3<> E_RRT;
	double mhalfInvSigmaSq = -0.5/(sigma*sigma);
	double cum_force_z = 0.;
	vector3<> cum_force(0.,0.,0.);
	for(Atom& atom: atoms)
		if(atom.sp == iSpecies)
		{
			//Compute wrapped Cartesian position:
			vector3<> posWrapped = atom.pos - center; //shift to origin of Gaussian
			for(int iDir=0; iDir<3; iDir++)
				posWrapped[iDir] -= floor(0.5 + posWrapped[iDir]); //to [-0.5, 0.5)
			vector3<> r = gInfo.R * posWrapped;
			
			//Extract appropriate reduced coordinate:
			double rSq = 0.; vector3<> rSq_r;
			rSq = r[2]*r[2]; 
			rSq_r = vector3<>(0., 0., 2*r[2]);
			
			//Accumulate energy, force and stress (if needed):
			double Ecur = U0 * exp(mhalfInvSigmaSq * rSq);
			vector3<> force = (-Ecur * mhalfInvSigmaSq) * rSq_r; //Cartesian force
			E += Ecur;
			cum_force += gInfo.RT * force; //convert to contravariant force
		}
	cum_force_z = cum_force[2];
	return cum_force_z;
}

double DynamicPlanarGaussianPotential::update_center(const GridInfo& gInfo, std::vector<Atom>& atoms, matrix3<>* E_RRTptr) const
{
	//Compute cumulative force in z-direction for all test centers, and update center to the test center closest to target force:
	//Note: this assumes the planar Gaussian is oriented normal to z-direction in Cartesian space
	

	double cum_force_z = cum_force_z(gInfo, atoms, E_RRTptr, &center);
	center[2] += step_size * cum_force_z;
	if(center[2] > 0.5) center[2] -= 1.0;
	if(center[2] < -0.5) center[2] += 1.0;
	return center[2];
}

