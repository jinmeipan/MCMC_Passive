Subroutine ObsModel(Theta,Tb)

Use CommonVars
Implicit None

Integer :: i
Real :: Profile(Nlyr,Npro)
Real,Intent(In) :: Theta(Ntheta)
Real,Intent(Out) :: Tb(Np,Nf)



Ctrl(1)=Nlyr

Xin(1)=Real(Nlyr)
Xin(2)=274.-Theta(Ntheta-2) !Soil temperature [K]
Xin(2)=Min(Xin(2),273.15)
Xin(3)=Theta(Ntheta-1)      !Soil water content [FRAC]
Xin(4)=Theta(Ntheta)        !Ground roughness [m]


Do i=1,Nlyr
  Profile(i,1)=Theta(iDz(i))  !Layer thickness [m]
  Profile(i,2)=Theta(iRho(i)) !Density [kg/m3]
  Profile(i,3)=Theta(iD(i))   !/1000.0   !Convert grain size to [m]
  Profile(i,4)=0.0            !Liquid fraction [m3/m3]
  Profile(i,5)=274.-Theta(iT(i))
End Do


! Call MEMLS, where Profile is snow properties, Tb_Bc in donward sky Tb,
! Xin contains the soil paramters, ScatOpt is the scattering coefficient option
! Opt_DmaxToPex is how to prepare exponential correlation length (directly use, for needs to be converted from snow grain size, Dmax)

Call SS_Model(Ctrl,Freq,Angle,Profile,Tb_Bc,Xin,Tb,ScatOpt,Opt_DmaxToPex)

End Subroutine ObsModel





