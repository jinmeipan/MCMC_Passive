%snowpit structure
classdef snowpit

properties
    %basic
    provider
    site
    year
    month
    date
    time
    %properties
	nlayer
    dz      %in m
    density %in kg/m^3
    T       %in K
    dmax    %in mm
    pex     %in mm
    mv      %m^3/m^3
    %soil properties
    mv_soil %m^3/m^3
    soilT   %K
    %tbs
    theta   %degree
    freq
    tbv     %K
    tbh     %K
    %statistics
    SD      %m
    SWE     %mm
    avg_density  %kg/m^3
    avg_T        %K
    avg_dmax     %mm
    avg_pex      %mm
    %note
    note
end

methods
    function snowpit=summary(snowpit)
        snowpit.SD=sum(snowpit.dz);
        snowpit.SWE=sum(snowpit.dz.*snowpit.density);
        weight=snowpit.dz.*snowpit.density;
        snowpit.avg_density=sum(weight.*snowpit.density)./sum(weight);
        snowpit.avg_T=sum(weight.*snowpit.T)./sum(weight);
        snowpit.avg_dmax=sum(weight.*snowpit.dmax)./sum(weight);
        snowpit.avg_pex=sum(weight.*snowpit.pex)./sum(weight);
    end
end

end