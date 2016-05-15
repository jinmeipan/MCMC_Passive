Subroutine ReadRunParams

Use CommonVars
Implicit None
Integer :: i,j
Real :: TMinLim,TMaxLim


Open(Unit=ParFileUnit,File=fnm(1),Status='Unknown')


! Number of pits to run (Npits)
Read(1,*)
Read(1,'(I5)') Npits
Print *, 'Npits=', Npits

! Number of variables in pit profile (Npro)
Read(1,*)
Read(1,'(I5)') Npro
Print *, 'Npro=', Npro

! Number of Theta variables (NthetaVars)
Read(1,*)
Read(1,'(I5)') NthetaVars
Print *, 'NthetaVars=', NthetaVars

! Number of iterations in the Markov Chain (Niter)
Read(1,*)
Read(1,'(I12)') Niter
Print *, 'Niter=', Niter

! Number of burn-in iterations in the Markov Chain (Njump)
Read(1,*)
Read(1,'(I6)') Nburn
Print *, 'Nburn=', Nburn

! Number of observation polarizations to use (Np)
Read(1,*)
Read(1,'(I5)') Np
Print *, 'Np=', Np

! Number of observation frequencies to use (Nf)
Read(1,*)
Read(1,'(I5)') Nf
Print *, 'Nf=', Nf

! Number of independent observation "times" (Nobs)
Read(1,*)
Read(1,'(I5)') Nobs
Print *, 'Nobs=', Nobs

! Number of Ctrl variables given to MEMLS (Nc)
Read(1,*)
Read(1,'(I5)') Nc

! Number of Auxiliary inputs given to MEMLS (Nx)
Read(1,*)
Read(1,'(I5)') Nx

! Error standard deviation of Tb observations (StdTb)
Read(1,*)
Read(1,'(F5.1)') StdTb
Print *, 'StdTb=', StdTb

! Scattering Coefficient Option 1 (Empirical), 2 (Born), or 3 (ScatOpt)
Read(1,*)
Read(1,'(I5)') ScatOpt
Print *, 'ScatOpt=', ScatOpt

! Method to convert from snow grain size to exponential correlation length Option 1 (Weissflujoch), 2 (Sodankyla), or 3 (AlreadyPex)
Read(1,*)
Read(1,'(I5)') Opt_DmaxToPex
Print *, 'Opt_DmaxToPex=', Opt_DmaxToPex

! Use prior information, 0 or 1. (UsePrior)
Read(1,*)
Read(1,'(I5)') UsePrior
Print *, 'UsePrior=', UsePrior

! Estimate dZ, -1, 0, 1  (EstimateDz)
Read(1,*)
Read(1,'(I5)') EstimateDz
Print *, 'EstimateDz=', EstimateDz

! Estimate rho, -1, 0, 1  (EstimateRho)
Read(1,*)
Read(1,'(I5)') EstimateRho
Print *, 'EstimateRho=', EstimateRho

! Estimate D, -1, 0, 1  (EstimateD)
Read(1,*)
Read(1,'(I5)') EstimateD
Print *, 'EstimateD=', EstimateD

! Estimate T, -1, 0, 1  (EstimateT)
Read(1,*)
Read(1,'(I5)') EstimateT
Print *, 'EstimateT=', EstimateT

! Estimate Soil Parameter (EstimateS)
Read(1,*)
Read(1,'(I5)') EstimateS
Print *, 'EstimateS=', EstimateS


! Observation frequencies (Freq), Nf lines
Allocate(Freq(Nf),Angle(Nf))
Read(1,*)
Do i=1,Nf
  Read(1,'(F5.1)') Freq(i)
End Do
Print *, 'Freq=', Freq

! Observation angles (Angle), Nf lines
Read(1,*)
Do i=1,Nf
  Read(1,'(F5.1)') Angle(i)
End Do
Print *, 'Angle=', Angle

! Tb boundary condition above snow surface (K), Nf*2 lines, V first, H second
Allocate(Tb_Bc(Np,Nf))
Read(1,*)
Do i=1,Np
  Do j=1,Nf
    Read(1,'(F11.1)') Tb_Bc(i,j)
  End Do
End Do
Print *, 'Tb_Bc=', Tb_Bc

! Minimum & maximum limits for layer thickness [m]
Read(1,*)
Read(1,'(F5.1)') DzMinLim
Read(1,'(F5.1)') DzMaxLim

! Minimum & maximum limits for density [kg/m3]
Read(1,*)
Read(1,'(F5.1)') RhoMinLim
Read(1,'(F5.1)') RhoMaxLim

! Minimum & maximum limits for grain diameter [mm]
Read(1,*)
Read(1,'(F5.1)') DMinLim
Read(1,'(F5.1)') DMaxLim


! Minimum & maximum limits for temperature [K]
Read(1,*)
Read(1,'(F6.1)') TMinLim
Read(1,'(F6.1)') TMaxLim

TpMaxLim=274-TMinLim !convert the temperature limit from K to 274-K used in the MCMC
TpMinLim=274-TMaxLim !convert the temperature limit from K to 274-K used in the MCMC


! Minimum & maximum limits for soil moisture [K]
Read(1,*)
Read(1,'(F6.1)') MvSMinLim
Read(1,'(F6.1)') MvSMaxLim

! Minimum & maximum limits for soil rms-height [m]
Read(1,*)
Read(1,'(F6.1)') GndSigMinLim
Read(1,'(F6.1)') GndSigMaxLim


Close(ParFileUnit)


End Subroutine ReadRunParams
