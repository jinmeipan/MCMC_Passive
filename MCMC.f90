Subroutine MCMC

Use CommonVars
Implicit None
Integer :: i,j,k
Integer :: Constrain_rho, Constrain_T

Real :: Z1(Niter*Nlyr),Z2(Niter*Nlyr),Z3(Niter*Nlyr),Z4(Niter*(Nlyr+1))
Real :: Z5(Niter),Z6(Niter)    !Z5 Z6 Added May21,2015
Real :: U1(Niter),U2(Niter),U3(Niter),U4(Niter)
Real :: U5(Niter), U6(Niter)   !U5 U6 Added May21,2015

Allocate(ThetaU(Ntheta),ThetaV(Ntheta))
Allocate(iDz(Nlyr),iRho(Nlyr),iD(Nlyr),iT(Nlyr+1))
Allocate(JmpDzStd(Nlyr),JmpRhoStd(Nlyr),JmpDStd(Nlyr),JmpTStd(Nlyr+1))
Allocate(iMvS(1),iGndSig(1),JmpMvSStd(1),JmpGndSigStd(1))


!Initialize, Set the first values of the state vector using the priors. Calculate the jump function before burn-in
Call Initialize


Call RandNormal(Niter*Nlyr,Z1,0+Pit+Nlyr+3) !Z*, To be multiplied by Jmp**Std to calculate jump step
Call RandNormal(Niter*Nlyr,Z2,100+Pit+Nlyr+5)
Call RandNormal(Niter*Nlyr,Z3,200+Pit+Nlyr+5)
Call RandNormal(Niter*(Nlyr+1),Z4,300+Pit+Nlyr+4)
Call RandUniform(Niter,U1,400+Pit+Nlyr+3)   !U*, used as be compared with likelihood ratio
Call RandUniform(Niter,U2,500+Pit+Nlyr+5)
Call RandUniform(Niter,U3,600+Pit+Nlyr+8)
Call RandUniform(Niter,U4,700+Pit+Nlyr+5)

!Added May21,2015
Call RandNormal(Niter,Z5,800+Pit+Nlyr+5)
Call RandNormal(Niter,Z6,900+Pit+Nlyr+5)
Call RandUniform(Niter,U5,1000+Pit+Nlyr+5)
Call RandUniform(Niter,U6,1100+Pit+Nlyr+5)


!Set if you want to constrain the density as smaller at surface, larger at bottom (increasing monotonically)
Constrain_rho=1

!Set if you want to constrain snow temperature as lower at surface, higher at bottom
Constrain_T=1





Do i=1,Niter

  !PRINT *, i,'/',Niter

  !Adjust the jump function when it reaches the Nburn
  If(i.Eq.Nburn+1) Call AdjustJump


  !Update thickness
  !If(EstimateDz.Eq.1) 
  Call Iterate(Z1((i-1)*Nlyr+1:i*Nlyr),U1(i),iDz,&
    JmpDzStd,DzMu,DzCov,Pu1,Na1,DzMinLim,DzMaxLim,Nlyr)


  !Update density
  If(Constrain_rho.eq.0)Then
        !If(EstimateRho.Eq.1) 
        Call Iterate(Z2((i-1)*Nlyr+1:i*Nlyr),U2(i),iRho,&
            JmpRhoStd,RhoMu,RhoCov,Pu2,Na2,RhoMinLim,RhoMaxLim,Nlyr)
  Else

        !Jinmei added, for this case, use the specific Iterate2 function
        !If(EstimateRho.Eq.1) Call Iterate_Density(Z2((i-1)*Nlyr+1:i*Nlyr),U2(i),iRho,&
        !    JmpRhoStd,RhoMu,RhoCov,Pu2,Na2,RhoMinLim,RhoMaxLim,Nlyr,i,U1(i)) !add U1(i), another any uniform distribution

        !If(EstimateRho.Eq.1) 
        Call Iterate2(Z2((i-1)*Nlyr+1:i*Nlyr),U2(i),iRho,&
            JmpRhoStd,RhoMu,RhoCov,Pu2,Na2,RhoMinLim,RhoMaxLim,Nlyr,i) !add U1(i), another any uniform distribution

  EndIf


  !Update grain size
  !If(EstimateD.Eq.1) 
    Call Iterate(Z3((i-1)*Nlyr+1:i*Nlyr),U3(i),iD,&
    JmpDStd,DMu,DCov,Pu3,Na3,DMinLim,DmaxLim,Nlyr)


  !Update temperature
  If(Constrain_T.eq.0)Then
        !If(EstimateT.Eq.1) 
        Call Iterate(Z4((i-1)*(Nlyr+1)+1:i*(Nlyr+1)),&
            U4(i),iT,JmpTStd,TMu,TCov,Pu4,Na4,TpMinLim,TpmaxLim,Nlyr+1)
  Else
        !Jinmei added, for this case, use the specific Iterate3 function
        !If(EstimateT.Eq.1) Call Iterate_T(Z4((i-1)*(Nlyr+1)+1:i*(Nlyr+1)),&
        !    U4(i),iT,JmpTStd,TMu,TCov,Pu4,Na4,TpMinLim,TpmaxLim,Nlyr+1,i,U3(i)) !add U3(i), another any uniform distribution

        !If(EstimateT.Eq.1) 
        Call Iterate4(Z4((i-1)*(Nlyr+1)+1:i*(Nlyr+1)),&
            U4(i),iT,JmpTStd,TMu,TCov,Pu4,Na4,TpMinLim,TpmaxLim,Nlyr+1,i) !add U3(i), another any uniform distribution
  EndIf

  !If(EstimateS.Eq.1)Then
    Call Iterate(Z5(i),U5(i),iMvS,JmpMvSStd,MvSMu,MvSCov,Pu5,Na5,&
        MvSMinLim,MvSMaxLim,1)
    Call Iterate(Z6(i),U6(i),iGndSig,JmpGndSigStd,GndSigMu,GndSigCov,&
       Pu6,Na6,GndSigMinLim,GndSigMaxLim,1)
  !Endif

  ThetaPost(1:Ntheta,i)=ThetaU

  Do j=1,Np
    Do k=1,Nf
      TbPost((j-1)*Nf+k,i)=TbU(j,k)
    End Do
  End Do

End Do

Deallocate(ThetaU,ThetaV,iDz,iRho,iD,iT,JmpDzstd,JmpRhoStd,JmpDStd,&
  JmpTStd,iMvS,iGndSig,JmpMvSStd,JmpGndSigStd)



End Subroutine MCMC




! Subroutine Iterate, to run the MCMC chain for general parameters in state vector
Subroutine Iterate(Z,U,iSel,JmpStd,Mu,Cov,Pu,Na,ThetaMin,ThetaMax,NiSel)

Use CommonVars
Implicit None

Integer,Intent(In) :: iSel(NiSel),NiSel
Integer,Intent(InOut) :: Na
Integer :: i
Real,Intent(In) :: Z(NiSel),JmpStd(NiSel),Mu(NiSel),Cov(NiSel),U,&
  ThetaMin,ThetaMax
Real,Intent(InOut) :: Pu
Real :: R,Pv

ThetaV=ThetaU


Do i=1,NiSel
  ThetaV(iSel(i))=ThetaU(iSel(i))+Z(i)*JmpStd(i)
  ThetaV(iSel(i))=Min(ThetaV(iSel(i)),ThetaMax)
  ThetaV(iSel(i))=Max(ThetaV(iSel(i)),ThetaMin)
End Do



Call ObsModel(ThetaV,TbV)

Call NProb(TbV,Fv)
Call LogNProb(ThetaV,iSel,Mu,Cov,Pv,NiSel)

If (UsePrior.Eq.0) Pv=Pu

R=Pv/Pu*Exp(-0.5*(Fv-Fu))

If (R.Gt.U) Then
  Na=Na+1
  ThetaU=ThetaV
  TbU=TbV
  Fu=Fv
  Pu=Pv
End If

End Subroutine Iterate




!!!!!!!! Jinmei added begin !!!!!!!!!
! Subroutine Iterate2, to run the MCMC chain for density, needs bottom density larger than surface density
Subroutine Iterate2(Z,U,iSel,JmpStd,Mu,Cov,Pu,Na,ThetaMin,ThetaMax,NiSel,Niter_temp)

Use CommonVars
Implicit None

Integer,Intent(In) :: iSel(NiSel),NiSel
Integer,Intent(InOut) :: Na
Integer :: i
Integer :: Niter_temp  !Nov25,2015 by Jinmei
Real,Intent(In) :: Z(NiSel),JmpStd(NiSel),Mu(NiSel),Cov(NiSel),U,&
ThetaMin,ThetaMax
Real,Intent(InOut) :: Pu
Real :: R,Pv
Real :: Temp_density(NiSel), density1, density2
Real :: Temp_Na2   !Nov25,2015 by Jinmei
Integer :: good_density

ThetaV=ThetaU

!added here:
Do i=1,NiSel
    Temp_density(i)=ThetaU(iSel(i))+Z(i)*JmpStd(i)
End Do


Temp_Na2=Real(Na)/Real(Niter_temp)



If(Temp_Na2.lt.0.20)Then
    density1=Temp_density(1)
    If (NiSel.Gt.1)then
        Do i=2,NiSel
            density2=Temp_density(i)
            If (density2>density1)then
                Temp_density(i)=density1
            Endif
            density1=Temp_density(i)
        Enddo
    Endif
    good_density=1

Else
    good_density=1
    If (NiSel.Gt.1)then
        Do i=2,NiSel
            density1=Temp_density(i-1)
            density2=Temp_density(i)
            If(density2>density1)Then
                good_density=0
            Endif
        Enddo
    Endif
Endif



Do i=1,NiSel
  ThetaV(iSel(i))=Temp_density(i)
  ThetaV(iSel(i))=Min(ThetaV(iSel(i)),ThetaMax)
  ThetaV(iSel(i))=Max(ThetaV(iSel(i)),ThetaMin)
End Do

Call ObsModel(ThetaV,TbV)

Call NProb(TbV,Fv)
Call LogNProb(ThetaV,iSel,Mu,Cov,Pv,NiSel)

If (UsePrior.Eq.0) Pv=Pu

R=Pv/Pu*Exp(-0.5*(Fv-Fu))


If (R.Gt.U .And. good_density.EQ.1) Then
    Na=Na+1
    ThetaU=ThetaV
    TbU=TbV
    Fu=Fv
    Pu=Pv
End If

End Subroutine Iterate2





Subroutine Iterate_Density(Z,U,iSel,JmpStd,Mu,Cov,Pu,Na,ThetaMin,ThetaMax,NiSel,Niter_temp,U2) !U2 added Dec22,2015

Use CommonVars
Implicit None

Integer,Intent(In) :: iSel(NiSel),NiSel
Integer,Intent(InOut) :: Na
Real,Intent(In) :: Z(NiSel),JmpStd(NiSel),Mu(NiSel),Cov(NiSel),U,ThetaMin,ThetaMax,U2
Real,Intent(InOut) :: Pu

Integer :: i
Integer :: Niter_temp  !Nov25,2015 by Jinmei
Real :: R,Pv
Real :: Temp_density(NiSel)
Real :: Z_Temp     !Dec21,2015


ThetaV=ThetaU

!added here:
Temp_density(1)=ThetaU(iSel(1))+Z(1)*JmpStd(1)

If (NiSel.Gt.1)Then
    Do i=2,NiSel
        Z_Temp=-abs(Z(i))/NiSel*0.5*U2*2.0 !Here use uniform distribution between 0 and 2..., and allow 0.5 as 60% percentile
        !Temp_density(i)=ThetaU(iSel(i))+Z_Temp*JmpStd(i)
        Temp_density(i)=Temp_density(i-1)+Z_Temp*JmpStd(i)
        !Print *,Z_Temp,Z(i),U2
    End Do
EndIf

Do i=1,NiSel
    ThetaV(iSel(i))=Temp_density(i)
    ThetaV(iSel(i))=Min(ThetaV(iSel(i)),ThetaMax)
    ThetaV(iSel(i))=Max(ThetaV(iSel(i)),ThetaMin)
End Do

Call ObsModel(ThetaV,TbV)

Call NProb(TbV,Fv)
Call LogNProb(ThetaV,iSel,Mu,Cov,Pv,NiSel)

If (UsePrior.Eq.0) Pv=Pu

R=Pv/Pu*Exp(-0.5*(Fv-Fu))


If (R.Gt.U) Then
    Na=Na+1
    ThetaU=ThetaV
    TbU=TbV
    Fu=Fv
    Pu=Pv
End If

End Subroutine Iterate_Density







! Subroutine Iterate4, to run the MCMC chain for 274-T, needs bottom 274-T smaller than surface 274-T
Subroutine Iterate4(Z,U,iSel,JmpStd,Mu,Cov,Pu,Na,ThetaMin,ThetaMax,NiSel,Niter_temp)

Use CommonVars
Implicit None

Integer,Intent(In) :: iSel(NiSel),NiSel
Integer,Intent(InOut) :: Na
Integer :: i
Integer :: Niter_temp  !Nov25,2015 by Jinmei
Real,Intent(In) :: Z(NiSel),JmpStd(NiSel),Mu(NiSel),Cov(NiSel),U,&
ThetaMin,ThetaMax
Real,Intent(InOut) :: Pu
Real :: R,Pv
Real :: Temp_Na4    !Nov25,2015 by Jinmei
Real :: Temp_T(NiSel), T1, T2
Integer :: good_T   !Nov25,2015 by Jinmei

ThetaV=ThetaU

!added here:
Do i=1,NiSel
    Temp_T(i)=ThetaU(iSel(i))+Z(i)*JmpStd(i)
End Do

Temp_Na4=Real(Na)/Real(Niter_temp)

If(Temp_Na4.lt.0.20)Then

    T1=Temp_T(1)
    If (NiSel.Gt.1)then
        Do i=2,NiSel-1  !Nov25,2015, Don't forget to exclude soil temperature from the constraining
            T2=Temp_T(i)
            If (T2<T1)then
                Temp_T(i)=T1
            Endif
            T1=Temp_T(i)
        Enddo
    Endif
    good_T=1

Else
    good_T=1
    If (NiSel.Gt.1)then
        Do i=2,NiSel-1  !Nov25,2015, Don't forget to exclude soil temperature from the constraining
            T1=Temp_T(i-1)
            T2=Temp_T(i)
            If(T2<T1)Then
                good_T=0
            Endif
        Enddo
    Endif
EndIf


Do i=1,NiSel
    ThetaV(iSel(i))=Temp_T(i)
    ThetaV(iSel(i))=Min(ThetaV(iSel(i)),ThetaMax)
    ThetaV(iSel(i))=Max(ThetaV(iSel(i)),ThetaMin)
End Do


Call ObsModel(ThetaV,TbV)

Call NProb(TbV,Fv)
Call LogNProb(ThetaV,iSel,Mu,Cov,Pv,NiSel)

If (UsePrior.Eq.0) Pv=Pu

R=Pv/Pu*Exp(-0.5*(Fv-Fu))


If(R.Gt.U .And. good_T.EQ.1)Then
Na=Na+1
ThetaU=ThetaV
TbU=TbV
Fu=Fv
Pu=Pv
End If

End Subroutine Iterate4





Subroutine Iterate_T(Z,U,iSel,JmpStd,Mu,Cov,Pu,Na,ThetaMin,ThetaMax,NiSel,Niter_temp,U2)  !U2 added Dec22,2015

Use CommonVars
Implicit None

Integer,Intent(In) :: iSel(NiSel),NiSel
Integer,Intent(InOut) :: Na
Real,Intent(In) :: Z(NiSel),JmpStd(NiSel),Mu(NiSel),Cov(NiSel),U,ThetaMin,ThetaMax,U2
Real,Intent(InOut) :: Pu

Integer :: i
Integer :: Niter_temp  !Nov25,2015 by Jinmei
Real :: R,Pv
Real :: Temp_T(NiSel)
Real :: Z_Temp     !Dec21,2015


ThetaV=ThetaU

!added here:
Temp_T(1)=ThetaU(iSel(1))+Z(1)*JmpStd(1)

If (NiSel.Gt.1)Then
    Do i=2,NiSel
        Z_Temp=abs(Z(i))/NiSel*0.5*U2*2.0  !Here use uniform distribution between 0 and 2..., and allow 0.5 as 60% percentile
        !Temp_T(i)=ThetaU(iSel(i))+Z_Temp*JmpStd(i)
        Temp_T(i)=Temp_T(i-1)+Z_Temp*JmpStd(i)
    End Do
EndIf


Do i=1,NiSel
    ThetaV(iSel(i))=Temp_T(i)
    ThetaV(iSel(i))=Min(ThetaV(iSel(i)),ThetaMax)
    ThetaV(iSel(i))=Max(ThetaV(iSel(i)),ThetaMin)
End Do


Call ObsModel(ThetaV,TbV)

Call NProb(TbV,Fv)
Call LogNProb(ThetaV,iSel,Mu,Cov,Pv,NiSel)

If (UsePrior.Eq.0) Pv=Pu

R=Pv/Pu*Exp(-0.5*(Fv-Fu))


If(R.Gt.U)Then
    Na=Na+1
    ThetaU=ThetaV
    TbU=TbV
    Fu=Fv
    Pu=Pv
End If

End Subroutine Iterate_T
!!!!!!!! Jinmei added finish !!!!!!!!!






! Subroutine Initialize, intialize the MCMC chain
Subroutine Initialize

Use CommonVars
Implicit None
Integer :: i
Integer :: AddAdjust

Na1=0
Na2=0
Na3=0
Na4=0
Na5=0 !Added May21,2015, needs to be added to CommonVars
Na6=0 !Added May21,2015, needs to be added to CommonVars

Do i=1,Nlyr
  iDz(i)=i
  iRho(i)=Nlyr+i
  iD(i)=2*Nlyr+i
End Do
Do i=1,Nlyr+1
  iT(i)=3*Nlyr+i
End Do

Do i=1,1
    iMvS(i)=4*Nlyr+2
    iGndSig(i)=4*Nlyr+3
End Do

Do i=1,Nlyr
  Select Case(EstimateDz)
    Case(-1)
      !ThetaU(iDz(i))=Exp(DzMu(i)+0.5*DzCov(i))
      ThetaU(iDz(i))=Exp(DzMu(i)-DzCov(i))
    Case(0)
      !ThetaU(iDz(i))=DzcmTrue(i)/100. !Convert to [m]
      ThetaU(iDz(i))=Exp(DzMu(i)+2.0*DzCov(i))
    Case(1)
      ThetaU(iDz(i))=Exp(DzMu(i)+0.5*DzCov(i))
  End Select


  Select Case(EstimateRho)
    Case(-1)
      !ThetaU(iRho(i))=Exp(RhoMu(i)+0.5*RhoCov(i))
      ThetaU(iRho(i))=Exp(RhoMu(i)-RhoCov(i))
    Case(0)
      !ThetaU(iRho(i))=RhoTrue(i)
      ThetaU(iRho(i))=Exp(RhoMu(i)+2.0*RhoCov(i))
    Case(1)
      ThetaU(iRho(i))=Exp(RhoMu(i)+0.5*RhoCov(i))
  End Select


  Select Case(EstimateD)
    Case(-1)
      !ThetaU(iD(i))=Exp(DMu(i)+0.5*DCov(i))
      ThetaU(iD(i))=Exp(DMu(i)-DCov(i))
    Case(0)
      !ThetaU(iD(i))=DmaxTrue(i)
      ThetaU(iD(i))=Exp(DMu(i)+2.0*DCov(i))
    Case(1)
      ThetaU(iD(i))=Exp(DMu(i)+0.5*DCov(i)) !这是正常的起点
  End Select
End Do



Select Case(EstimateT)
  Case(-1)
    Do i=1,Nlyr+1
      !ThetaU(iT(i))=Exp(TMu(i)+0.5*TCov(i))
      ThetaU(iT(i))=Exp(TMu(i)-TCov(i))
    End Do
  Case(0)
    Do i=1,Nlyr+1
      !ThetaU(iT(i))=274.-TsnowTrue(i)
      ThetaU(iT(i))=Exp(TMu(i)+2.0*TCov(i))
    End Do
    !ThetaU(iT(Nlyr+1))=274.-TgTrue
  Case(1)
    Do i=1,Nlyr+1
      ThetaU(iT(i))=Exp(TMu(i)+0.5*TCov(i))
    End Do
End Select



Select Case(EstimateS)
  Case(-1)
    !ThetaU(iMvS(1))=Exp(MvSMu(1)+0.5*MvSCov(1))
    !ThetaU(Ntheta)=Exp(GndSigMu(1)+0.5*GndSigCov(1))
    ThetaU(iMvS(1))=Exp(MvSMu(1)-MvSCov(1))
    ThetaU(Ntheta)=Exp(GndSigMu(1)-GndSigCov(1))
  Case(0)
    !ThetaU(iMvS(1))=MvSTrue
    !ThetaU(Ntheta)=GndSigTrue
    ThetaU(iMvS(1))=Exp(MvSMu(1)+2.0*MvSCov(1))
    ThetaU(Ntheta)=Exp(GndSigMu(1)+2.0*GndSigCov(1))
  Case(1)
    ThetaU(iMvS(1))=Exp(MvSMu(1)+0.5*MvSCov(1))
    ThetaU(Ntheta)=Exp(GndSigMu(1)+0.5*GndSigCov(1))
End Select


AddAdjust=1

Do i=1,Nlyr
  JmpDzStd(i)=Sqrt( (Exp(DzCov(i))-1.)*Exp(2*DzMu(i)+DzCov(i)) )&
    /Sqrt(Real(Nlyr*Ntheta+2))/10./Sqrt(Real(Nobs)) *1.8
    !/changed from Nlyr*3 to Nlyr*Nthta+2, and from 10. to 5., Jul1
  JmpRhoStd(i)=Sqrt( (Exp(RhoCov(i))-1.)*Exp(2*RhoMu(i)+RhoCov(i)) )&
    /Sqrt(Real(Nlyr*Ntheta+2))/10./Sqrt(Real(Nobs)) *2.0
  JmpDStd(i)=Sqrt( (Exp(DCov(i))-1.)*Exp(2*DMu(i)+DCov(i)) )&
    /Sqrt(Real(Nlyr*Ntheta+2))/10./Sqrt(Real(Nobs)) *2.8


  !Jinmei add, to distinguish when using good prior or bad prior
  If(AddAdjust.eq.1)Then
        !to limit the Jump functon from too small
        If(JmpDzStd(i)<0.02*1.8)Then
            JmpDzStd(i)=Sqrt( (Exp(DzCov(i))-1.)*Exp(2*DzMu(i)+DzCov(i)) )&
                /Sqrt(Real(Nlyr*Ntheta+2))/5./Sqrt(Real(Nobs)) *1.8
        Endif

        If(JmpRhoStd(i)<3*2.0)Then
            JmpRhoStd(i)=Sqrt( (Exp(RhoCov(i))-1.)*Exp(2*RhoMu(i)+RhoCov(i)) )&
                /Sqrt(Real(Nlyr*Ntheta+2))/5./Sqrt(Real(Nobs)) *2.0
        Endif

        If(JmpDStd(i)<0.001*2.8)Then
            JmpDStd(i)=Sqrt( (Exp(DCov(i))-1.)*Exp(2*DMu(i)+DCov(i)) )&
                /Sqrt(Real(Nlyr*Ntheta+2))/5./Sqrt(Real(Nobs)) *2.8
        Endif

        !to limit the jump function from too big, Jul2, Jinmei
        If(JmpDzStd(i)>0.05*1.8)Then  !Avg jump of thickness no Larger than 5 cm
            JmpDzStd(i)=0.05*1.8
        Endif

        If(JmpRhoStd(i)>20.0*2.0)Then !Avg jump of density no larger than 20.0 kg/m^3
            JmpRhoStd(i)=20.0*2.0
        Endif

        If(JmpDStd(i)>0.004*2.8)Then  !Avg jump of dmax no Larger than 0.004 mm
            JmpDStd(i)=0.004*2.8
        Endif
   EndIf
End Do


Do i=1,Nlyr+1
  JmpTStd(i)=Sqrt( (Exp(TCov(i))-1.)*Exp(2*TMu(i)+TCov(i)) )&
    /Sqrt(Real(Nlyr*Ntheta+2))/10./Sqrt(Real(Nobs)) *6.0

   !Jinmei added, to distinguish when using good prior or bad prior
  If(AddAdjust.eq.1)Then
        !to limit the Jump functon from too small
        If(JmpTStd(i)<0.3*6.0)Then
            JmpTStd(i)=Sqrt( (Exp(TCov(i))-1.)*Exp(2*TMu(i)+TCov(i)) )&
                /Sqrt(Real(Nlyr*Ntheta+2))/5./Sqrt(Real(Nobs)) *6.0
        Endif

        !to limit the jump function from too big, Jul2, Jinmei
        If(JmpTStd(i)>0.5*6.0)Then  !Avg jump of T no larger than 0.5 K
            JmpTStd(i)=0.5 *6.0
        Endif
   Endif
End Do



JmpMvSStd(1)=Sqrt( (Exp(MvSCov(1))-1.)*Exp(2*MvSMu(1)+MvSCov(1)) )&
/Sqrt(Real(Nlyr*Ntheta+2))/10./Sqrt(Real(Nobs)) *16.0

JmpGndSigStd(1)=Sqrt( (Exp(GndSigCov(1))-1.)*Exp(2*GndSigMu(1)+GndSigCov(1)) )&
/Sqrt(Real(Nlyr*Ntheta+2))/10./Sqrt(Real(Nobs)) *16.0


! Jinmei add, revised according to unit change
  If(AddAdjust.eq.1)Then
        ! to prevent the jump function from too big
        If(JmpMvSStd(1)>0.003*16.0)Then     !Avg jump of soil moisture no larger than 0.3%
            JmpMvSStd(1)=0.003*16.0
        Endif

        If(JmpGndSigStd(1)>0.0001*16.0)Then !Avg jump of ground roughness no larger than 0.1 mm
            JmpGndSigStd(1)=0.0001*16.0
        Endif
 EndIf


Print *,'Initialized jump dz (cm):',JmpDzStd*100
Print *,'Initialized jump rho (kg/m^3):',JmpRhoStd
Print *,'Initialized jump D (mm):',JmpDStd
Print *,'Initialized jump T (C):',JmpTStd
Print *,'Initialized jump MvS (%):',JmpMvSStd*100.0
Print *,'Initialized jump GndSig (mm):',JmpGndSigStd*1000.0


Call ObsModel(ThetaU,TbU)

Call NProb(TbU,Fu)
Call LogNProb(ThetaU,iDz,DzMu,DzCov,Pu1,Nlyr)
Call LogNProb(ThetaU,iRho,RhoMu,RhoCov,Pu2,Nlyr)
Call LogNProb(ThetaU,iD,DMu,DCov,Pu3,Nlyr)
Call LogNProb(ThetaU,iT,TMu,TCov,Pu4,Nlyr+1)
Call LogNProb(ThetaU,iMvS,MvSMu,MvSCov,Pu5,1)
Call LogNProb(ThetaU,iGndSig,GndSigMu,GndSigCov,Pu6,1)

End Subroutine Initialize









Subroutine AdjustJump

Use CommonVars
Implicit None
Integer :: i
Integer :: AddAdjust
Real :: DzStd(Nlyr),RhoStd(Nlyr),DStd(Nlyr),TStd(Nlyr+1)
Real :: MvSStd,GndSigStd !May21,2015 added



AddAdjust=1


Do i=1,Nlyr
    Call CalcStdDev(ThetaPost(iDz(i),1:Nburn),Nburn,DzStd(i))
    JmpDzStd(i)=DzStd(i)/Sqrt(Real(3*Nlyr))*2.38/2.0 *1.8

    Call CalcStdDev(ThetaPost(iRho(i),1:Nburn),Nburn,RhoStd(i))
    JmpRhoStd(i)=RhoStd(i)/Sqrt(Real(3*Nlyr))*2.38 *2.0

    Call CalcStdDev(ThetaPost(iD(i),1:Nburn),Nburn,DStd(i))
    JmpDStd(i)=DStd(i)/Sqrt(Real(3*Nlyr))*2.38 *2.8


    If(AddAdjust.eq.1)Then
        If(JmpDzStd(i)>0.05*1.8)Then  !Avg jump of thickness no larger than 5 cm
            JmpDzStd(i)=0.05*1.8
        Endif

        If(JmpRhoStd(i)>20.0*2.0)Then !Avg jump of density no larger than 20 kg/m3
            JmpRhoStd(i)=20.0*2.0
        Endif

        If(JmpDStd(i)>0.004*2.8)Then  !Avg jump of grain size no larger than 0.004 mm
            JmpDStd(i)=0.004*2.8
        Endif
    Endif
End Do



Do i=1,Nlyr+1
    Call CalcStdDev(ThetaPost(iT(i),1:Nburn),Nburn,TStd(i))
    JmpTStd(i)=TStd(i)/Sqrt(Real(3*Nlyr))*2.38 *6.0


    If(AddAdjust.eq.1)Then
        If(JmpTStd(i)>0.5*6.0)Then  !Avg jump of temperature no larger than 0.5 K
            JmpTStd(i)=0.5*6.0
        Endif
    Endif
End Do


!Added May21,2015
Call CalcStdDev(ThetaPost(Ntheta-1,1:Nburn),Nburn,MvSStd)
JmpMvSStd=MvSStd/Sqrt(Real(3*Nlyr))*2.38/2.0 *16.0

Call CalcStdDev(ThetaPost(Ntheta,1:Nburn),Nburn,GndSigStd)
JmpGndSigStd=GndSigStd/Sqrt(Real(3*Nlyr))*2.38/2.0 * 16.0


If(AddAdjust.eq.1)Then
    If(JmpMvSStd(1)>0.003*16.0)Then !Avg jump of soil moisture no larger than 0.3%
        JmpMvSStd(1)=0.003*16.0
    Endif

    If(JmpGndSigStd(1)>0.0001*16.0)Then  !Avg jump of soil roughness no larger than 0.1 mm
        JmpGndSigStd(1)=0.0001*16.0
    Endif
Endif

Print *,'Adjuested jump dz (cm):',JmpDzStd*100
Print *,'Adjuested jump rho (kg/m^3):',JmpRhoStd
Print *,'Adjuested jump D (mm):',JmpDStd
Print *,'Adjuested jump T (C):',JmpTStd
Print *,'Adjuested jump MvS (%):',JmpMvSStd*100.0
Print *,'Adjuested jump GndSig (mm):',JmpGndSigStd*1000.0

End Subroutine AdjustJump













Subroutine LogNProb(Theta,iSel,Mu,Cov,P,NiSel)

Use CommonVars
Implicit None
Integer,Intent(In) :: iSel(NiSel),NiSel
Integer :: i
Real,Intent(In) :: Theta(Ntheta),Mu(Nlyr),Cov(Nlyr)
Real,Intent(Out) :: P

P=0.0

!Do i=1,Nlyr !Jinmei revised May21,2015, needs to check
Do i=1,NiSel
  P=P+(Log(Theta(iSel(i)))-Mu(i))**2./Cov(i)
End Do
P=Exp(-0.5*P)

!Do i=1,Nlyr
Do i=1,NiSel
  P=P/Theta(iSel(i))
End Do

End Subroutine LogNProb





Subroutine NProb(Tb,F)

Use CommonVars
Implicit None
Integer :: i,j,k
Real,Intent(In) :: Tb(Np,Nf)
Real,Intent(Out) :: F

F=0.
Do i=1,Np
  Do j=1,Nf
    Do k=1,Nobs
      F=F+(TbObs(i,j,k)-Tb(i,j))**2/StdTb**2.
    End Do
  End Do
End Do

End Subroutine NProb






Subroutine CalcStdDev(X,N,StdDevX)

Implicit None
Integer,Intent(In) :: N
Integer :: i
Real,Intent(In) :: X(N)
Real,Intent(Out) :: StdDevX
Real :: MuX

MuX=Sum(X)/Real(N)
StdDevX=0.
Do i=1,N
  StdDevX=StdDevX+(X(i)-MuX)**2
End Do
StdDevX=Sqrt(StdDevX/Real(N-1))

End Subroutine CalcStdDev
