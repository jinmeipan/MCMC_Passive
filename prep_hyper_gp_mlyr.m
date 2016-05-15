%this document using the different priors for different layrs

%write hyper.txt
%
% ---  Discussion on how to determin lognormal mean and std -----
% m and v are mean and varianc of the data X
%m = 1;
%v = 2;
% mu and sigma are mean and standard reviation of log(X)
%mu = log((m^2)/sqrt(v+m^2));
%sigma = sqrt(log(v/(m^2)+1));
%[exp(mu-sigma),exp(mu),exp(mu+sigma)]=
%0.202411359667921       0.577350269189626           1.6468113937884;
%[m-sqrt(v),m,m+sqrt(v)]=
%-0.414213562373095                         1          2.41421356237309;
%if the first run starts from exp(mu+sigma), than this method is
%acceptable
% ------------------------------------------------------------------


%inputs:
%-gpr: in the format of structure prior2  !!a different strucutre from prior!!
%-model: determine to input pex or dmax prior


%density, grain size, temperature, thickness
function prep_hyper_gp_mlyr(site,model,month,folder0,file)

switch site
    case 'sdkl'
        load('/Users/office/Desktop/CLPX_data_FromNSIDC/prior_sdkl_morelyr.mat');
        switch month
            case 11
                gpr=prior_sdkl_nov;
            case 12
                gpr=prior_sdkl_dec;
            case 1
                gpr=prior_sdkl_jan;
            case 2
                gpr=prior_sdkl_feb;
            case 3
                gpr=prior_sdkl_mar;
            case 4
                gpr=prior_sdkl_apr;
        end
        
    case 'churchill'
        load('/Users/office/Desktop/CLPX_data_FromNSIDC/prior_churchill_morelyr.mat');
        switch month
            case 1
                gpr=prior_churchill_jan;
            case 2
                gpr=prior_churchill_feb;
            case 3
                gpr=prior_churchill_mar;
            case 4
                gpr=prior_churchill_apr;
        end
        
    case 'clpx'
            load('/Users/office/Desktop/CLPX_data_FromNSIDC/prior_clpx_morelyr.mat');
            gpr=prior_clpx_feb;
end


%soil moisture
MvS_mu=-2.872302235;
MvS_sigma=0.832554611;

%soil roughness
GndSig_mu=-4.951743776;   
GndSig_sigma=0.832554611;



%% 
fid_hyper=fopen([folder0,'temp/hyperpar_',file,'.txt'],'w');
Max_lyr=6;
for ilyr=1:Max_lyr
    
    fprintf(fid_hyper,'%2i\n',ilyr);

    
    %thickness
    fprintf(fid_hyper, 'Layer thickness (log-normal),m \n');
    
    swe_mean=gpr.swe_mean;
    
    percent=0.4;
    swe_std=swe_mean*percent;

    density_mean=gpr.density_mean{1}(1);
    density_std=gpr.density_std{1}(1);
    
    [sd_mean,sd_std]=swe2sd(swe_mean,swe_std,density_mean,density_std);
    
    M=sd_mean;
    V=sd_std.^2;
    
    M_plan=M/ilyr;
	V_plan=V/ilyr;
    dz_mu = log((M_plan^2)/sqrt(V_plan+M_plan^2));
    dz_sigma = sqrt(log(V_plan/(M_plan^2)+1));    
        
    for k=1:ilyr
        fprintf(fid_hyper, '%9.6f\n',dz_mu);
    end
    for k=1:ilyr
        fprintf(fid_hyper, '%9.6f\n',dz_sigma.^2);
    end
      
    
%     fprintf(fid_hyper, 'Layer thickness (log-normal),m \n');
%     
%     for k=1:ilyr
%         fprintf(fid_hyper, '%9.6f\n',gpr.dz_mu{ilyr}(k));
%     end
%     for k=1:ilyr
%         fprintf(fid_hyper, '%9.6f\n',gpr.dz_sigma{ilyr}(k).^2);
%     end
    
    
    %grain size
    fprintf(fid_hyper, 'Grain size (log-normal),mm \n');
    
    for k=1:ilyr
        switch model
            case 'hut'
                fprintf(fid_hyper, '%9.6f\n',gpr.dmax_mu{ilyr}(k));
            otherwise
                fprintf(fid_hyper, '%9.6f\n',gpr.pex_mu{ilyr}(k));
        end
    end
    for k=1:ilyr
        switch model
            case 'hut'
                fprintf(fid_hyper, '%9.6f\n',gpr.dmax_sigma{ilyr}(k).^2);
            otherwise
                fprintf(fid_hyper, '%9.6f\n',gpr.pex_sigma{ilyr}(k).^2);
        end
    end
    
    
    %density
    fprintf(fid_hyper, 'Density (log-normal),kg/m^3 \n');
    
    for k=1:ilyr
        fprintf(fid_hyper, '%9.6f\n',gpr.density_mu{ilyr}(k));
    end
    for k=1:ilyr
        fprintf(fid_hyper, '%9.6f\n',gpr.density_sigma{ilyr}(k).^2);
    end
    
    
    %temperature
    fprintf(fid_hyper, '274.-Temperature (log-normal),K \n');
    
    for k=1:ilyr
        fprintf(fid_hyper, '%9.6f\n',gpr.T_mu{ilyr}(k)); %snow T from bottom to surface
    end
    fprintf(fid_hyper, '%9.6f\n',gpr.T_mu{ilyr}(end));	 %ground T
    
    
    for k=1:ilyr
        fprintf(fid_hyper, '%9.6f\n',gpr.T_sigma{ilyr}(k).^2);
    end
    fprintf(fid_hyper, '%9.6f\n',gpr.T_sigma{ilyr}(end).^2);
    
    
    %soil moisture and roughness, revised for R2
    fprintf(fid_hyper, 'Soil volumetric water content,FRAC \n');
    fprintf(fid_hyper, '%9.6f\n',MvS_mu);       %ground T
    fprintf(fid_hyper, '%9.6f\n',MvS_sigma.^2);
    
    
	fprintf(fid_hyper, 'Soil surface rougness,m \n');
    fprintf(fid_hyper, '%9.6f\n',GndSig_mu);       %ground T
    fprintf(fid_hyper, '%9.6f\n',GndSig_sigma.^2);
end

fclose(fid_hyper);

end