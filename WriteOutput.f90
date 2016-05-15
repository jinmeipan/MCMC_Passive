Subroutine WriteOutput

Use CommonVars
Implicit None
Integer :: i,j


Do i=1,Ntheta
  Write(ThetaOutFileUnit) (ThetaPost(i,j),j=1,Niter)
  !Print *,'(ThetaPost(i,j)',ThetaPost(i-1,j-1)
End Do

Do i=1,Np*Nf
  Write(TbOutFileUnit) (TbPost(i,j),j=1,Niter)
  !Print *,'(TbPost(i,j)',TbPost(i-1,j-1)
End Do

Write(AcptFileUnit,'(I5,1X,I5)') Pit,Nlyr
Write(AcptFileUnit,'(F5.3)') Real(Na1)/Real(Niter)
Write(AcptFileUnit,'(F5.3)') Real(Na2)/Real(Niter)
Write(AcptFileUnit,'(F5.3)') Real(Na3)/Real(Niter)
Write(AcptFileUnit,'(F5.3)') Real(Na4)/Real(Niter)
Write(AcptFileUnit,'(F5.3)') Real(Na5)/Real(Niter)
Write(AcptFileUnit,'(F5.3)') Real(Na6)/Real(Niter)


End Subroutine WriteOutput
