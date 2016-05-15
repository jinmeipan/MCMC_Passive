Subroutine ReadObs

Use CommonVars
Implicit None
Integer :: i,j,k

Read(ObsFileUnit,*) 

! Read TB Observations from tb_obs.txt file. The file contains a line of data
! First Nf ones are for vertical polarization, Second Nf ones are for horizontal 
! polarization

Do i=1,Np
  Do j=1,Nf
    Do k=1,Nobs
      Read(ObsFileUnit,'(F5.1)') TbObs(i,j,k)
    End Do
  End Do
End Do
End Subroutine ReadObs
