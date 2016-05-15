%Run params prepare
%input:
%-model: which observation model; e.g. model='hut'; where 'emp'=empirical memls; 'iba'='iba','hut'=hut model combination
%-site: where is the pits; e.g. site='sdkl'
%-sp: the snowpit information in type of Structure Snowpit
%-stdTb: observation error of TB (K); eg. stdTb=2.0


function [freq_index,freq,theta,TskyStr]=prep_runparam2(sp,model,site,npol,stdTb,folder0,file,freq_plan)

npits=size(sp,1);

switch model
    %memls options
    case 'emp'
     scatopt=1; %use empirical memls
    case 'iba'
     scatopt=2; %use iba
    case 'hut'
     scatopt=2; %use combination of Halli & Roy
    otherwise
        return
        disp('please choose method to calculate ks');
end
disp(['scatopt=',num2str(scatopt)]);



%
load('/Users/office/Google Drive Present/!Solve_HUT&MEMLS/Tsky.mat')
%includes parameter TskyALL, with
%rows: 6.925, 10.65, 18.7, 36.5 and 90 GHz
%columns:from left to right, 50, 53 and 55 deg.
%Christ's data is at 53degree,6.9, 19, 37 and 89 GHz
% TskyAll=TskyAll*0.0+2.7; %comment


switch site
    case 'sdkl'
        if(strcmp(freq_plan,'4freq')==1)
            freq_index=[1,2,3,4];
            freq=[10.65,18.7,36.5,89];
            TskyStr=TskyAll([2,3,4,5],1);
        end
        if(strcmp(freq_plan,'3freq')==1)
            freq_index=[1,3,4];
            freq=[10.65,36.5,89];
            TskyStr=TskyAll([2,4,5],1);
        end
        if(strcmp(freq_plan,'2freq')==1)
            freq_index=[2,3];
            freq=[18.7,36.5];
            TskyStr=TskyAll([3,4],1);
        end
        
        theta=50;  %50 for sodankyla, 53 for churchill, 55 for clpx
    case 'churchill'
        if(strcmp(freq_plan,'4freq')==1)
            freq_index=[1,2,3,4];
            freq=[6.925,18.7,36.5,89];
            TskyStr=TskyAll([1,3,4,5],2);
        end
        if(strcmp(freq_plan,'3freq')==1)
            freq_index=[1,3,4];
            freq=[6.925,36.5,89];
            TskyStr=TskyAll([1,4,5],2);
        end
        if(strcmp(freq_plan,'2freq')==1)
            freq_index=[2,3];
            freq=[18.7,36.5];
            TskyStr=TskyAll([3,4],2);
        end
        
        theta=53;
    case 'clpx'
        freq_index=[1,2,3];
        freq=[18.7,36.5,89];
        theta=55;
        TskyStr=TskyAll([3,4,5],3);

%         freq_index=[1,2,3];
%         freq=[10.65,18.7,36.5,89];
%         theta=55;
%         TskyStr=TskyAll([2,3,4,5],1);
        if(strcmp(freq_plan,'2freq')==1)
            freq_index=[1,2];
            freq=[18.7,36.5];
            theta=55;
            TskyStr=TskyAll([3,4],3);
        end
        
end
disp(['theta=',num2str(theta)])
disp(['freq=',num2str(freq)])

%Method to convert from dmax to pex: 1 (Weissflujoch), 2 (Sodankyla), or 3 (Directly use for both HUT and pex)\n');
pexopt=3;  %to be modified (actually fixed already)


%%
% npol=1;
% stdTb=2.0;


%% begin to write
Niter=10000;  %default values=10000
Niter=Niter*3;  %to be revised!
Estimate=[1,1,1,1,1];     %to be revised!
% Estimate=[0,0,0,0,0];     %to be revised!
% Estimate=[-1,-1,-1,-1,-1];     %to be revised!
% Estimate=[0, 1, 1,1,1];     %to be revised!
% Estimate=[-1,1,1,1,1];     %to be revised!


fid_rp=fopen([folder0,'temp/RunParamsV2_',file,'.txt'],'w');

fprintf(fid_rp,'Number of pits to run (Npits)\n');
fprintf(fid_rp,'%3i\n',npits);
fprintf(fid_rp,'Number of variables in pit profile (Npro) (Don?t change)\n');
fprintf(fid_rp,'5\n');
fprintf(fid_rp,'Number of Theta variables (NthetaVars)\n');
fprintf(fid_rp,'4\n');
fprintf(fid_rp,'Number of iterations in the Markov Chain (Niter)\n');
%fprintf(fid_rp,'10000\n');
fprintf(fid_rp,'%8i\n',Niter);
fprintf(fid_rp,'Number of burn-in iterations in the Markov Chain (Nburn)\n');
fprintf(fid_rp,'2000\n');
fprintf(fid_rp,'Number of observation polarizations to use (Np)\n');
fprintf(fid_rp,'%3i\n',npol);
fprintf(fid_rp,'Number of observation frequencies to use (Nf)\n');
fprintf(fid_rp,'%3i\n',length(freq));
fprintf(fid_rp,'Number of independent observation "times" (Nobs)\n');
fprintf(fid_rp,'1\n');
fprintf(fid_rp,'Number of Ctrl variables given to MEMLS (Nc), 4 (Don?t change)\n');
fprintf(fid_rp,'4\n');
fprintf(fid_rp,'Number of Auxiliary inputs given to MEMLS (Nx), 8 (Don?t change)\n');
fprintf(fid_rp,'8\n');
fprintf(fid_rp,'!Error standard deviation of Tb observations (StdTb)\n');
fprintf(fid_rp,'%4.1f\n',stdTb);
fprintf(fid_rp,'Scattering Coefficient Option (ScatOpt): 1 (Empirical, Halli), 2 (Born, HUT combination), or 3 (Roy)\n');
fprintf(fid_rp,'%3i\n',scatopt);
fprintf(fid_rp,'Method to convert from snow grain size to exponential correlation length Option (Opt_DmaxToPex): 1 (Weissflujoch), 2 (Sodankyla), or 3 (AlreadyPex)\n');
fprintf(fid_rp,'%3i\n',pexopt);
fprintf(fid_rp,'Use prior information, 0 or 1. (UsePrior)\n');
fprintf(fid_rp,'1\n');
fprintf(fid_rp,'Estimate dZ, -1, 0, 1  (EstimateDz)\n');
fprintf(fid_rp,'%3i\n',Estimate(1));
fprintf(fid_rp,'Estimate rho, -1, 0, 1  (EstimateRho)\n');
fprintf(fid_rp,'%3i\n',Estimate(2));
fprintf(fid_rp,'Estimate D, -1, 0, 1  (EstimateD)\n');
fprintf(fid_rp,'%3i\n',Estimate(3));
fprintf(fid_rp,'Estimate T, -1, 0, 1  (EstimateT)\n');
fprintf(fid_rp,'%3i\n',Estimate(4));
fprintf(fid_rp,'Estimate S, -1, 0, 1  (EstimateS)\n');
fprintf(fid_rp,'%3i\n',Estimate(5));
fprintf(fid_rp,'Observation frequencies (Freq), Nf lines\n');
for i=1:length(freq)
	fprintf(fid_rp,'%6.2f\n',freq(i));
end
fprintf(fid_rp,'Observation angles (Angle), Nf lines\n');
for i=1:length(freq)
	fprintf(fid_rp,'%6.2f\n',theta);
end
fprintf(fid_rp,'Tb boundary condition above snow surface (K), Nf*2 lines, V first, H second\n');
%for i=1:length(freq)*npol
%    fprintf(fid_rp,'2.7\n');
%end
for j=1:npol
    for i=1:length(freq)
        fprintf(fid_rp,'%8.3f\n',TskyStr(i));
    end
end
fprintf(fid_rp,'Minimum & maximum limits for layer thickness [m]\n');
fprintf(fid_rp,'0.01\n'); %revised to prevent errors!
fprintf(fid_rp,'10.0\n');
fprintf(fid_rp,'Minimum & maximum limits for density [kg/m3]\n');
fprintf(fid_rp,'50.0\n');
fprintf(fid_rp,'917.0\n');
fprintf(fid_rp,'Minimum & maximum limits for grain diameter [mm]\n');
fprintf(fid_rp,'0.02\n'); %revised from 0.1mm, because this is also limit for pex
fprintf(fid_rp,'5.0\n');
fprintf(fid_rp,'Minimum & maximum limits for temperature [K]\n');
fprintf(fid_rp,'243.15\n');
fprintf(fid_rp,'273.15\n');
fprintf(fid_rp,'Minimum & maximum limits for soil moisture [Frac]\n');
fprintf(fid_rp,'0\n');
fprintf(fid_rp,'1\n');
fprintf(fid_rp,'Minimum & maximum limits for soil rms-height [m]\n');
fprintf(fid_rp,'0\n');
fprintf(fid_rp,'0.1\n');
fclose(fid_rp);

clear pexopt i fid_rp j


end