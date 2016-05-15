%input:
%-sp: in the format of snowpit structure


function prep_truesp(sp,model,folder0)

for i=1:length(sp);
    sd_c(i,1)=sp(i).SD;
    swe_c(i,1)=sp(i).SWE;
    density_c(i,1)=sp(i).avg_density;
    T_c(i,1)=sp(i).avg_T;
    Tsoil_c(i,1)=sp(i).T(1);
    switch model
        case 'hut'
            gs_c(i,1)=sp(i).avg_dmax;
        otherwise
            gs_c(i,1)=sp(i).avg_pex;
    end
end

output=[sd_c;nan;nan;swe_c;nan;nan;density_c;nan;nan;gs_c;nan;nan;274-T_c;nan;nan;274-Tsoil_c];
% open output


clear i sd_c swe_c density_c
clear T_c Tsoil_c gs_c

save([folder0,'temp/sp.mat'],'sp')


end

%appendix A:
% properties of snowpits structure
%     %basic
%     provider
%     site
%     year
%     month
%     date
%     time
%     %properties
% 	nlayer
%     dz      %in m
%     density %in kg/m^3
%     T       %in K
%     dmax    %in mm
%     pex     %in mm
%     mv      %m^3/m^3
%     %soil properties
%     mv_soil %m^3/m^3
%     soilT   %K
%     %tbs
%     theta   %degree
%     freq
%     tbv     %K
%     tbh     %K
%     %statistics
%     SD      %m
%     SWE     %mm
%     avg_density  %kg/m^3
%     avg_T        %K
%     avg_dmax     %mm
%     avg_pex      %mm
%     %note
%     note
% end