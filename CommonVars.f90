Module CommonVars

Implicit None



! A. Index of input-output files
Integer,Parameter :: ParFileUnit=1,ObsFileUnit=2,TrueFileUnit=3,&
  HyperFileUnit=4,ThetaOutFileUnit=7,AcptFileUnit=8,TbOutFileUnit=9


! B. Basic parameters from RunParams.txt
Integer :: Npits,Npro,NthetaVars,Niter,Nburn,Np,Nf,Nobs,Nc,Nx
Integer :: UsePrior     !Use Prior or not to calculate likelihood function
Integer :: EstimateDz,EstimateRho,EstimateT,EstimateD,EstimateS !Exclude some parameters in MCMC estimation if setted as zero
Integer :: ScatOpt
Integer :: Opt_DmaxToPex

Real :: StdTb !Standard deviation of Tb observations
Real :: DzMinLim,DzMaxLim,RhoMinLim,RhoMaxLim,DMinLim,DMaxLim
Real :: TpMinLim,TpMaxLim !274-TMaxLim,274-TMinLim
Real :: MvSMinLim,MvSMaxLim,GndSigMinLim,GndSigMaxLim !Min and Max limit of params

Real,Dimension(:),Allocatable :: Freq,Angle
Real,Dimension(:,:,:),Allocatable :: TbObs !TB observations
Real,Dimension(:,:),Allocatable :: Tb_Bc   !Downward Tb at snow surface


Real,Dimension(:),Allocatable :: DzMu,DzCov,RhoMu,RhoCov,DMu,DCov,TMu,TCov  !Medium and Std of snow parameter priors
Real,Dimension(:),Allocatable :: MvSMu,MvSCov,GndSigMu,GndSigCov !Medium and Std of soil parameter priors


! C. Intermediate Paramters
Integer :: Pit  !The index of current snowpit
Integer :: Ntheta !Total number of parameters to be iterated (do MCMC run) for one snowpit
Integer :: Nlyr !Number of layers in snowpack, to be changed for different layer plans
Integer :: Na1,Na2,Na3,Na4 !Number of acceptance
Integer :: Na5,Na6 !Added May21,2015

Integer,Dimension(:),Allocatable :: iDz,iRho,iD,iT,iMvS,iGndSig !Indexes of snow paramers in Theta, Added May21,2015
Integer,Dimension(:),Allocatable :: Ctrl    !Array to save Dimensions
Real,Dimension(:),Allocatable :: Xin        !Array to save soil parameters


! Ctrl vector
! Ctrl(1): Number of snow layers, Nlyr
! Ctrl(2): Number of polarizations, Np
! Ctrl(3): Number of snow parameters, Npro
! Ctrl(4): Number of frequencies,Nf

! Xin vector speficies 8 auxliliary inputs for MEMLS, mainly related to soil
! Xin(1): Number of snow layers, Nlyr
! Xin(2): Soil Temperature, Theta(Ntheta-2)
! Xin(3): Soil volumbetric water content [FRAC], Theta(Ntheta-1)
! Xin(4): Ground roughness [m], Theta(Ntheta)
! Xin(5): Soil density [kg/m^3], fixed values for MEMLS run but could be revised in SetupMEMLS function
! Xin(6): Soil sand content [%],fixed values for MEMLS run but could be revised in SetupMEMLS function
! Xin(7): Soil silt content [%],fixed values for MEMLS run but could be revised in SetupMEMLS function
! Xin(8): Soil clay content [%],fixed values for MEMLS run but could be revised in SetupMEMLS function



Real :: Fu,Fv
Real :: Pu1,Pu2,Pu3,Pu4,Pu5,Pu6
!Real :: Pv1,Pv2,Pv3,Pv4,Pv5,Pv6
        !Fu,the posterior probablity of simulated Tb
        !Pu*,the posterior probablity of snow parameters
        !Pv is not use actually

Real,Dimension(:),Allocatable :: ThetaU,ThetaV
        !Vectors of theta in MCMC run, where U is the state vector for the past iteration, 
        !V is the state vector for the current iteration (after likelyhood function judgement)
Real,Dimension(:),Allocatable :: JmpDzStd,JmpRhoStd,JmpDStd,JmpTStd !Jump function
Real,Dimension(:),Allocatable :: JmpMvSStd,JmpGndSigStd
Real,Dimension(:,:),Allocatable :: TbU,TbV
!Real,Dimension(:,:),Allocatable :: Tbs


! D. Output arrays
Real,Dimension(:,:),Allocatable :: ThetaPost,TbPost



! E. Default values of the two newly-added soil parameters
Real :: MvSTrue, GndSigTrue



! F. Added May24,2015, used to save filenames
Integer :: iFile
Character*150 fnm(7)



! The following variables are not used Now!!
!Integer :: iLyr
!Real :: TgTrue     !Default values of other snow parameters
!Real,Dimension(:),Allocatable :: Fliq,Sppt,Pex
!Real,Dimension(:),Allocatable :: DzcmTrue,RhoTrue,TsnowTrue,DmaxTrue
!Real,Dimension(:),Allocatable :: Tsnow,Rho,Dzcm,Dmax
!Real,Dimension(:,:),Allocatable :: TrueProfile



End Module



