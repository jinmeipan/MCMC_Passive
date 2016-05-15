Subroutine ReadHyperPar
! HyperPar revised to include two soil parameters, soil moisture [FRAC] and soil roughness rms-height [m]

Use CommonVars
Implicit None
Integer :: i


Read(HyperFileUnit,*) !Skip for number of layers


! Read layer thickness prior, from snow bottom to snow surface
Read(HyperFileUnit,*)
Do i=1,Nlyr
  Read(HyperFileUnit,'(F9.6)') DzMu(i)
End Do

Do i=1,Nlyr
  Read(HyperFileUnit,'(F9.6)') DzCov(i)
End Do


! Read snow grain size prior, from snow bottom to snow surface
Read(HyperFileUnit,*)
Do i=1,Nlyr
  Read(HyperFileUnit,'(F9.6)') DMu(i)
End Do
Do i=1,Nlyr
  Read(HyperFileUnit,'(F9.6)') DCov(i)
End Do


! Read snow density prior, from snow bottom to snow surface
Read(HyperFileUnit,*)
Do i=1,Nlyr
  Read(HyperFileUnit,'(F9.6)') RhoMu(i)
End Do
Do i=1,Nlyr
  Read(HyperFileUnit,'(F9.6)') RhoCov(i)
End Do


! Read snow temperature prior, from snow bottom to snow surface
! and the last one is the ground temperature
Read(HyperFileUnit,*)
Do i=1,Nlyr+1
  Read(HyperFileUnit,'(F9.6)') TMu(i)
End Do
Do i=1,Nlyr+1
  Read(HyperFileUnit,'(F9.6)') TCov(i)
End Do

! Read soil moisture prior
Read(HyperFileUnit,*)
Do i=1,1
    Read(HyperFileUnit,'(F9.6)') MvSMu(i)
    Read(HyperFileUnit,'(F9.6)') MvSCov(i)
End Do


! Read soil roughness rms-height prior
Read(HyperFileUnit,*)
Do i=1,1
    Read(HyperFileUnit,'(F9.6)') GndSigMu(i)
    Read(HyperFileUnit,'(F9.6)') GndSigCov(i)
End Do



End Subroutine ReadHyperPar
