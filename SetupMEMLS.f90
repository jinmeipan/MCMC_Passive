Subroutine SetupMEMLS

Use CommonVars
Implicit None

Allocate(Ctrl(Nc),Xin(Nx)) !Never deallocated

!Static MEMLS Inputs
Ctrl(2)=Np
Ctrl(3)=Npro
Ctrl(4)=Nf

!Fixed soil parameters
Xin(5)=1500.0  !Soil bulk density in kg/m3
Xin(6)=30.63 !Sand content of soil texture [%]
Xin(7)=55.89 !Silt content of soil texture [%]
Xin(8)=13.48 !Clay content of soil texture [%]

MvSMaxLim=1.0-Xin(5)/2650.0 !change the maximum volumetric water content from 1 to porosity


!Set this if you do not want to estimate the soil moisture and roughness
If (EstimateS==0) Then
  MvSTrue=0.08      !Soil volumetric water content [FRAC], Set soil moisture as 8%
  GndSigTrue=0.001  !Soil roughness [m], Set defaul roughness as 1 mm
Endif

End Subroutine SetupMEMLS
