Program MetroMEMLS

Use CommonVars

Implicit None
Integer :: CountI,CountF

Character*150 foo
Character*150 mypath
Character*150 tempfile
Integer :: ipath,j


Call System_Clock(CountI)

Print *, "MetroMEMLS version 2.0, 21 Dec 2015"


!Read filenames
If(.true.)Then
    Call getarg(0,foo)
    ipath=index(foo,'/',.true.)
    mypath=foo(1:ipath)
    Print *, mypath
    Open(81,file=trim(mypath)//'FILENAME.txt',status='old')
Else
    Open(81,file='/Users/office/Downloads/MCMC/FILENAME.txt',status='old')
EndIf


Do iFile=1,7
    Read(81,'(a)') fnm(iFile)
    If(.true.)Then
        tempfile=trim(mypath)//trim(fnm(iFile))
        Print *,tempfile
        fnm(iFile)=tempfile
    Endif
End Do
Close(81)


Open(Unit=ObsFileUnit,File=fnm(2),Status='Old')
!Open(Unit=TrueFileUnit,File=fnm(4),Status='Old')
!Read(TrueFileUnit,*)
!Read(TrueFileUnit,*)
Open(Unit=ThetaOutFileUnit,File=fnm(7),Status='Unknown',Form='Unformatted')
Open(Unit=TbOutFileUnit,File=fnm(6),Status='Unknown',Form='Unformatted')
Open(Unit=AcptFileUnit,File=fnm(5),Status='Unknown')


Call ReadRunParams

Call SetupMEMLS

Allocate(TbObs(Np,Nf,Nobs),TbU(Np,Nf),TbV(Np,Nf))


Do Pit=1,Npits

  Print *, 'Running pit #',Pit,'/',Npits

  Call ReadObs
  !Read(TrueFileUnit,*)
  !Open(Unit=HyperFileUnit,File='!hyperpar.txt',Status='Old')

  Open(Unit=HyperFileUnit,File=fnm(3),Status='Old')


  Do Nlyr=1,2  !Revised, run for 1 to 6 layers only

    Print *, 'Layer plan #',Nlyr

    Ntheta=Nlyr*NthetaVars+3  !Revised 18/5/15, add two more params: soil moisture and roughness

    !Allocate(TrueProfile(Nlyr,Npro))
    !Allocate(DzcmTrue(Nlyr),RhoTrue(Nlyr),DmaxTrue(Nlyr),TsnowTrue(Nlyr))
    Allocate(DzMu(Nlyr),DzCov(Nlyr),RhoMu(Nlyr),RhoCov(Nlyr),&
      DMu(Nlyr),DCov(Nlyr),TMu(Nlyr+1),TCov(Nlyr+1))
    Allocate(MvSMu(1),MvSCov(1),GndSigMu(1),GndSigCov(1))
    Allocate(ThetaPost(Ntheta,Niter),TbPost(Np*Nf,Niter))

    !Call SetupTrueLayerData

    Call ReadHyperPar

    Call MCMC

    !Print *,(ThetaPost(1,j),j=1,Niter)
    !Print *,(TbPost(1,j),j=1,Niter)

    Call WriteOutput


    Deallocate(DzMu,DzCov,RhoMu,RhoCov,DMu,DCov,TMu,TCov)
    Deallocate(MvSMu,MvSCov,GndSigMu,GndSigCov)
    Deallocate(ThetaPost,TbPost)

    !Deallocate(TrueProfile)
    !Deallocate(Sppt,Fliq,Pex)
    !Deallocate(DzcmTrue,RhoTrue,DmaxTrue,TsnowTrue)

  End Do
  Close(HyperFileUnit)
End Do


Close(ObsFileUnit)
!Close(TrueFileUnit)
Close(ThetaOutFileUnit)
Close(AcptFileUnit)


Deallocate(Freq,Angle,TbObs,Tb_Bc)
Deallocate(Ctrl,Xin)
Deallocate(TbU,TbV)


Call System_Clock(CountF)
Print *, 'Job took',Real(CountF-CountI)/1000./60.,' minutes.'


End Program MetroMEMLS


