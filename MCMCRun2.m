classdef MCMCRun2
    
properties

%run parameters
UsePrior
%sizes
Npits
NumThetaVars 
Niter
Nburn
Np
Nf
Nobs
%Tb observationss
TbObs
%hyperparameters
StdTb
DzMu
DzCov
RhoMu
RhoCov
DMu
DCov
TMu
TCov
MvSMu
MvSCov
GndSigMu
GndSigCov
%The general experiment output
theta_post
acceptance
TbPost
%chains for all parameters: chains{Nlyr,Npits}(1:Nlyr,1:Niter)
chains_dz
chains_rho
chains_D
chains_T
chains_Tsoil
chains_mvs
chains_sig
chains_sd
chains_swe
%mean of the chains for each layer after Burn-in: mean{Nlyr,Npits}(1:Niter)
%this mean is not the directly taking the mean of the chains, but taking
%the exp of the mean of the log of the mean.
md_dz
md_rho
md_D
md_T
md_Tsoil
md_mvs
md_sig
md_sd
md_swe
%standard derivation of the chains after Burn-in: std{Nlyr,Npits}(1:Niter)
%this standard derivation is not directly taking the std, but the std error
%around the mean above.
std_dz
std_rho
std_D
std_T
std_Tsoil
std_mvs
std_sig
std_sd
std_swe
%the final output !!!!
J %intermediate variable for model selection
nHat %number of snow layers
sdHat
SWEHat
rhoavgHat    
DavgHat
TavgHat
TsoilHat
MvSHat
SigHat
sdStdHat
SWEStdHat
%true pit data  !!!!
nTrue    %number of layers of true snowpits
sdTrue    %snow depth of true snowpits
SWETrue    %get SWE of true snowpits
TruePits   %snowpit structure
%filenames
filename
filename_runparam
filename_tbobs
filename_hyperpar
filename_truetheta
filename_acceptance
filename_posttb
filename_posttheta
%User settings
Max_Lyr     %maximum number of layers
Obs_Model   %observation model type
SturmClass  %snow type
lyrplan     %layer plan for the final output
nlyr_choose %number of layers if we fix the values to pick
run_location %if on osc super computer (=1), or no local desktop(=0)
%Limits, Jinmei Dec26,2015
DzMinLim
DzMaxLim
RhoMinLim
RhoMaxLim
DMinLim
DMaxLim
TMinLim
TMaxLim
MvSMinLim
MvSMaxLim
GndSigMinLim
GndSigMaxLim
end %properties




methods    
    %1) constructor
    %2) reads
    %3) analysis
    %4) utilities
    %5) plotting
    %6) writing
    
function Run=Main(Run)
   
    
   %Run.filename='/Users/office/Downloads/MCMC/MCMCRunData_V2/FILENAME_T.txt';
   length=Run.filename;
   if(length==0)
       error('please show me where the filename is at first!');
   end
   
   %Set basic running options
%    Run.Max_Lyr=6;
   
   %Layer plans
%    Run.lyrplan=1; %1=use modelSelection to select number of layers; 2=contrained number of layers to 2
%    Run.nlyr_choose=2; %write when lyplan=2
   

   %% Read output
   Run=Run.ReadFilenames;
   Run=Run.ReadRunParams;
   Run=Run.ReadObs;
   Run=Run.ReadHyperPar;
   Run=Run.ReadAcceptance; 
   
 
   Run=Run.ReadTbPost;
   Run=Run.ReadThetaOut;
   Run=Run.CalcThetaOut;
   
   
   %Select the best result
   Run=Run.ModelSelection;
   
   
   %Plot output
   %Run=Run.PlotAcceptance(Run.MaxLyr);
   %Run=Run.PlotThetaOut(Run,Pit,MaxLyr);
   %Run=Run.PlotThetaOutHist(Run,Pit,MaxLyr);
   %Run.PlotThetaOutHist(3,2)
   %Run=Run.PlotThetaCompare(Run,Pit);
end




%% read inputs & outputs
function Run=ReadFilenames(Run)
    fid=fopen(Run.filename);
    Run.filename_runparam=fscanf(fid,'%s \n',1);
    Run.filename_tbobs=fscanf(fid,'%s \n',1);
    Run.filename_hyperpar=fscanf(fid,'%s \n',1);
    Run.filename_truetheta=fscanf(fid,'%s \n',1);
    Run.filename_acceptance=fscanf(fid,'%s \n',1);
    Run.filename_posttb=fscanf(fid,'%s \n',1);
    Run.filename_posttheta=fscanf(fid,'%s \n',1);
    fclose(fid);
    
    %return complete names of subfiles;
    %we know that there are all under the same folder of FILENAME.txt
	string=Run.filename_tbobs;
    [pathstr,name,ext]=fileparts(Run.filename);
    Run.filename_runparam=[pathstr,'/',Run.filename_runparam];
    Run.filename_tbobs=[pathstr,'/',Run.filename_tbobs];
    Run.filename_hyperpar=[pathstr,'/',Run.filename_hyperpar];
    Run.filename_truetheta=[pathstr,'/',Run.filename_truetheta];
    Run.filename_acceptance=[pathstr,'/',Run.filename_acceptance];
    Run.filename_posttb=[pathstr,'/',Run.filename_posttb];
    Run.filename_posttheta=[pathstr,'/',Run.filename_posttheta];
    
    
    
    %Read truesnowpits
    string=Run.filename_tbobs;
    idx=strfind(string,'tb_obs_');
    string2=string(idx+7:end-4); %example of string2: 'sdmaremp', or 'sdmar'
    
    temp2=length(string2);
    if(temp2>7)
        complexname=1;
        
        %add the model names
        string3=string2(end-2:end); %example of string3, 'emp'
%         Run.Obs_Model=1; %1=MEMLS,2=HUT, to be revised!
        string2=string2(1:end-3);
        
        switch string3
            case 'iba'
                Run.Obs_Model=1;
            case 'emp'
                Run.Obs_Model=1;
            case 'hut'
                Run.Obs_Model=2;
            otherwise
                Run.Obs_Model=1;
                disp(['the observation model info is not given,', ... 
                    'set it as MEMLS; then pex will be extrated from truepits intead of max']);
        end
    end
    
    
    %set snow class
    Run.SturmClass=1; %1='Taiga',2='Alpine', to be revised!
    
    switch string2
        case 'sdnov'
            load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_sdkl.mat');
            Run.TruePits=tsp_sdkl_nov;
            clear tsp_sdkl_nov tsp_sdkl_dec tsp_sdkl_jan tsp_sdkl_feb tsp_sdkl_mar
        case 'sddec'
            load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_sdkl.mat');
            Run.TruePits=tsp_sdkl_dec;
            clear tsp_sdkl_nov tsp_sdkl_dec tsp_sdkl_jan tsp_sdkl_feb tsp_sdkl_mar
        case 'sdjan'
            load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_sdkl.mat');
            Run.TruePits=tsp_sdkl_jan;
            clear tsp_sdkl_nov tsp_sdkl_dec tsp_sdkl_jan tsp_sdkl_feb tsp_sdkl_mar
        case 'sdfeb'
            load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_sdkl.mat');
            Run.TruePits=tsp_sdkl_feb;
            clear tsp_sdkl_nov tsp_sdkl_dec tsp_sdkl_jan tsp_sdkl_feb tsp_sdkl_mar
        case 'sdmar'
            load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_sdkl.mat');
            Run.TruePits=tsp_sdkl_mar;
            clear tsp_sdkl_nov tsp_sdkl_dec tsp_sdkl_jan tsp_sdkl_feb tsp_sdkl_mar
        case 'chjan'
            load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_churchill.mat');
            Run.TruePits=tsp_churchill_jan;
            clear tsp_churchill_jan tsp_churchill_feb tsp_churchill_mar tsp_churchill_apr
        case 'chfeb'
            load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_churchill.mat');
            Run.TruePits=tsp_churchill_feb;
            clear tsp_churchill_jan tsp_churchill_feb tsp_churchill_mar tsp_churchill_apr
        case 'chmar'
            load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_churchill.mat');
            Run.TruePits=tsp_churchill_mar;
            clear tsp_churchill_jan tsp_churchill_feb tsp_churchill_mar tsp_churchill_apr
        case 'chapr'
            load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_churchill.mat');
            Run.TruePits=tsp_churchill_apr;
            clear tsp_churchill_jan tsp_churchill_feb tsp_churchill_mar tsp_churchill_apr
        case 'clpxfeb'
            load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_clpx_lsos.mat');
            Run.TruePits=sp_clpx_lsos;
            clear sp_clpx_lsos
            Run.SturmClass=2; %1='Tundra' or 'Taiga',2='Alpine', to be revised!
    end
    
end %function ReadFilenames


function Run=ReadRunParams(Run)
   fid=fopen(Run.filename_runparam,'r');
   fgetl(fid); Run.Npits=fscanf(fid,'%f \n',1);
   fgetl(fid); fgetl(fid);            
   fgetl(fid); Run.NumThetaVars=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.Niter=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.Nburn=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.Np=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.Nf=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.Nobs=fscanf(fid,'%f \n',1);
   fgetl(fid); fgetl(fid); 
   fgetl(fid); fgetl(fid); 
   fgetl(fid); Run.StdTb=fscanf(fid,'%f \n',1);
   fgetl(fid); fgetl(fid);
   fgetl(fid); fgetl(fid);
   fgetl(fid); Run.UsePrior=fscanf(fid,'%f \n',1);
   %Estimate
   fgetl(fid); fgetl(fid);
   fgetl(fid); fgetl(fid);
   fgetl(fid); fgetl(fid);
   fgetl(fid); fgetl(fid);
   fgetl(fid); fgetl(fid);
   %Observed.freq
   fgetl(fid); 
   for i=1:Run.Nf
       fgetl(fid);
   end
   %angle
   fgetl(fid); 
   for i=1:Run.Nf
       fgetl(fid);
   end
   %tsky
   fgetl(fid); 
   for i=1:Run.Nf*Run.Np
       fgetl(fid);
   end
   %limits
   fgetl(fid);Run.DzMinLim=fscanf(fid,'%f \n',1); Run.DzMaxLim=fscanf(fid,'%f \n',1);
   fgetl(fid);Run.RhoMinLim=fscanf(fid,'%f \n',1); Run.RhoMaxLim=fscanf(fid,'%f \n',1);
   fgetl(fid);Run.DMinLim=fscanf(fid,'%f \n',1); Run.DMaxLim=fscanf(fid,'%f \n',1);
   fgetl(fid);Run.TMinLim=fscanf(fid,'%f \n',1); Run.TMaxLim=fscanf(fid,'%f \n',1);
   fgetl(fid);Run.MvSMinLim=fscanf(fid,'%f \n',1); Run.MvSMaxLim=fscanf(fid,'%f \n',1);
   fgetl(fid);Run.GndSigMinLim=fscanf(fid,'%f \n',1); Run.GndSigMaxLim=fscanf(fid,'%f \n',1);
   fclose(fid);
   
   temp1=Run.TMinLim;
   temp2=Run.TMaxLim;
   Run.TMinLim=274-temp2;
   Run.TMaxLim=274-temp1;
   clear temp1 temp2
end %function ReadRunParams


function Run=ReadObs(Run)
    %read in the observations
    Nchan=Run.Np*Run.Nf;
    Run.TbObs=nan(Run.Nobs,Nchan,Run.Npits);
    fid=fopen(Run.filename_tbobs,'r');
    for i=1:Run.Npits,        
        fgetl(fid);
        for j=1:Nchan,
            Run.TbObs(:,j,i)=fscanf(fid,'%f \n',Run.Nobs)';
        end
    end
    fclose(fid);
end %function ReadObs


function Run=ReadHyperPar(Run)
    MaxNlyr=Run.Max_Lyr;
    fid=fopen(Run.filename_hyperpar,'r');
    for Nlyr=1:MaxNlyr,
        fgetl(fid);
        fgetl(fid);       
        Run.DzMu{Nlyr}=fscanf(fid,'%f \n',Nlyr);
        Run.DzCov{Nlyr}=fscanf(fid,'%f \n',Nlyr);       
        fgetl(fid);
        Run.DMu{Nlyr}=fscanf(fid,'%f \n',Nlyr);
        Run.DCov{Nlyr}=fscanf(fid,'%f \n',Nlyr);
        fgetl(fid);
        Run.RhoMu{Nlyr}=fscanf(fid,'%f \n',Nlyr);
        Run.RhoCov{Nlyr}=fscanf(fid,'%f \n',Nlyr);        
        fgetl(fid);
        Run.TMu{Nlyr}=fscanf(fid,'%f \n',Nlyr+1);
        Run.TCov{Nlyr}=fscanf(fid,'%f \n',Nlyr+1);        
        fgetl(fid);
        Run.MvSMu{Nlyr}=fscanf(fid,'%f \n',1);
        Run.MvSCov{Nlyr}=fscanf(fid,'%f \n',1);
	    fgetl(fid);
        Run.GndSigMu{Nlyr}=fscanf(fid,'%f \n',1);
        Run.GndSigCov{Nlyr}=fscanf(fid,'%f \n',1);
    end
    fclose(fid);
end


function Run=ReadAcceptance(Run)
   MaxNlyr=Run.Max_Lyr;
   Run.acceptance.dz=nan(Run.Npits,MaxNlyr);
   Run.acceptance.rho=nan(Run.Npits,MaxNlyr);
   Run.acceptance.D=nan(Run.Npits,MaxNlyr);
   Run.acceptance.T=nan(Run.Npits,MaxNlyr);
   Run.acceptance.MvS=nan(Run.Npits,MaxNlyr);
   Run.acceptance.GndSig=nan(Run.Npits,MaxNlyr);
   
   fid=fopen(Run.filename_acceptance,'r');

   for i=1:Run.Npits,
       for j=1:MaxNlyr,
           fgetl(fid);           
           Run.acceptance.dz(i,j)=fscanf(fid,'%f \n',1);
           Run.acceptance.rho(i,j)=fscanf(fid,'%f \n',1);
           Run.acceptance.D(i,j)=fscanf(fid,'%f \n',1);
           Run.acceptance.T(i,j)=fscanf(fid,'%f \n',1);
           Run.acceptance.MvS(i,j)=fscanf(fid,'%f \n',1);
           Run.acceptance.GndSig(i,j)=fscanf(fid,'%f \n',1);
       end
   end
   fclose(fid);
end %function ReadAcceptance



function Run=ReadTbPost(Run)    
    NlyrMax=Run.Max_Lyr;
    Nchan=Run.Np*Run.Nf;
    Run.TbPost=nan(Nchan,Run.Niter,Run.Npits,NlyrMax);
    
    %read in the Tb output
    fid=fopen(Run.filename_posttb,'rb');

    for i=1:Run.Npits,    
        for Nlyr=1:NlyrMax,
            for j=1:Run.Np*Run.Nf,
                fseek(fid,4,'cof');
                if(Run.run_location==1)
                    fseek(fid,4,'cof');
                end
                A=fread(fid,Run.Niter,'float');
                Run.TbPost(j,1:Run.Niter,i,Nlyr)=A;
                fseek(fid,4,'cof');
                if(Run.run_location==1)
                    fseek(fid,4,'cof');
                end
%                 open A
            end
        end      
    end
    fclose(fid);    
end %function ReadTbPost



function Run=ReadThetaOut(Run)
    for Nlyr=1:Run.Max_Lyr,
        Ntheta=Nlyr*Run.NumThetaVars+3;  %now it is right, to be revised
        Run.theta_post{Nlyr}=nan(Ntheta,Run.Niter,Run.Npits);
    end

    fid=fopen(Run.filename_posttheta,'rb');

    for i=1:Run.Npits,    
        for Nlyr=1:Run.Max_Lyr, 
            Ntheta=Nlyr*Run.NumThetaVars+3; %now it is right, to be revised
            for j=1:Ntheta, 
                fseek(fid,4,'cof');
                
                if(Run.run_location==1)
                    fseek(fid,4,'cof');
                end 
%                 Run.theta_post{Nlyr}(j,1:Run.Niter,i)=...
%                     fread(fid,Run.Niter,'float'); 
                
                A=fread(fid,Run.Niter,'float');
                Run.theta_post{Nlyr}(j,1:Run.Niter,i)=A;
                fseek(fid,4,'cof');
                if(Run.run_location==1)
                    fseek(fid,4,'cof');
                end
            end
        end      
    end
    
    
    fclose(fid);
end%function ReadThetaOut



%% try to deal with the output...
function [iDz,iRho,iD,iT,iTsoil,iMvS,iGndSig] = GetPostIndices(Run,n)    
    iDz=1:n;
    iRho=n+1:n*2;
    iD=n*2+1:n*3;
    iT=n*3+1:n*4;
    iTsoil=n*4+1;
    iMvS=n*4+2;
    iGndSig=n*4+3;
end %subfunction GetPostIndices


function Run=CalcThetaOut(Run)
    for i=1:Run.Npits
        %summarize chains
        for Nlyr=1:Run.Max_Lyr
            [iDz,iRho,iD,iT,iTsoil,iMvS,iGndSig] = Run.GetPostIndices(Nlyr);
            
            Run.chains_dz{Nlyr,i}=Run.theta_post{Nlyr}(iDz,:,i); %Nlyr rows, Niter columns
            Run.chains_rho{Nlyr,i}=Run.theta_post{Nlyr}(iRho,:,i);
            Run.chains_D{Nlyr,i}=Run.theta_post{Nlyr}(iD,:,i);
            Run.chains_T{Nlyr,i}=Run.theta_post{Nlyr}(iT,:,i);
            Run.chains_Tsoil{Nlyr,i}=Run.theta_post{Nlyr}(iTsoil,:,i);
            Run.chains_mvs{Nlyr,i}=Run.theta_post{Nlyr}(iMvS,:,i);
            Run.chains_sig{Nlyr,i}=Run.theta_post{Nlyr}(iGndSig,:,i);
            
            if(Nlyr>1)
                Run.chains_sd{Nlyr,i}=sum(Run.chains_dz{Nlyr,i},1);
                Run.chains_swe{Nlyr,i}=sum(Run.chains_dz{Nlyr,i}.*Run.chains_rho{Nlyr,i},1);
            else
                Run.chains_sd{Nlyr,i}=Run.chains_dz{Nlyr,i};
                Run.chains_swe{Nlyr,i}=Run.chains_dz{Nlyr,i}.*Run.chains_rho{Nlyr,i};
            end
            
            %summarize mean and std of the chains after burn-in (!!the
            %very original is mode, the most frequent values)
            %here revised it again by using the assumptions of log-normal
            %distributions
            iBurn=Run.Nburn+1:Run.Niter;
%             iBurn=Run.Nburn+1:20000; %to be revised
            
%             Run.md_dz{Nlyr,i}=mean(Run.chains_dz{Nlyr,i}(:,iBurn),2);
%             Run.std_dz{Nlyr,i}=std(Run.chains_dz{Nlyr,i}(:,iBurn),0,2);
            [logmean,logstd]=calc_lognormal_mean(Run.chains_dz{Nlyr,i}(:,iBurn),[Run.DzMinLim,Run.DzMaxLim]);
            Run.md_dz{Nlyr,i}=logmean;
            Run.std_dz{Nlyr,i}=logstd;
            
%             Run.md_rho{Nlyr,i}=mean(Run.chains_rho{Nlyr,i}(:,iBurn),2);
%             Run.std_rho{Nlyr,i}=std(Run.chains_rho{Nlyr,i}(:,iBurn),0,2);
            [logmean,logstd]=calc_lognormal_mean(Run.chains_rho{Nlyr,i}(:,iBurn),[Run.RhoMinLim,Run.RhoMaxLim]);
            Run.md_rho{Nlyr,i}=logmean;
            Run.std_rho{Nlyr,i}=logstd;


%             Run.md_D{Nlyr,i}=mean(Run.chains_D{Nlyr,i}(:,iBurn),2);
%             Run.std_D{Nlyr,i}=std(Run.chains_D{Nlyr,i}(:,iBurn),0,2);
            [logmean,logstd]=calc_lognormal_mean(Run.chains_D{Nlyr,i}(:,iBurn),[Run.DMinLim,Run.DMaxLim]);
            Run.md_D{Nlyr,i}=logmean;
            Run.std_D{Nlyr,i}=logstd;


%             Run.md_T{Nlyr,i}=mean(Run.chains_T{Nlyr,i}(:,iBurn),2);
%             Run.std_T{Nlyr,i}=std(Run.chains_T{Nlyr,i}(:,iBurn),0,2);
            [logmean,logstd]=calc_lognormal_mean(Run.chains_T{Nlyr,i}(:,iBurn),[Run.TMinLim,Run.TMaxLim]);
            Run.md_T{Nlyr,i}=logmean;
            Run.std_T{Nlyr,i}=logstd;


%             Run.md_Tsoil{Nlyr,i}=mean(Run.chains_Tsoil{Nlyr,i}(:,iBurn),2);
%             Run.std_Tsoil{Nlyr,i}=std(Run.chains_Tsoil{Nlyr,i}(:,iBurn),0,2);
            [logmean,logstd]=calc_lognormal_mean(Run.chains_Tsoil{Nlyr,i}(:,iBurn),[Run.TMinLim,Run.TMaxLim]);
            Run.md_Tsoil{Nlyr,i}=logmean;
            Run.std_Tsoil{Nlyr,i}=logstd;
            
            
%             Run.md_mvs{Nlyr,i}=mean(Run.chains_mvs{Nlyr,i}(:,iBurn),2);
%             Run.std_mvs{Nlyr,i}=std(Run.chains_mvs{Nlyr,i}(:,iBurn),0,2);
            [logmean,logstd]=calc_lognormal_mean(Run.chains_mvs{Nlyr,i}(:,iBurn),[Run.MvSMinLim,Run.MvSMaxLim]);
            Run.md_mvs{Nlyr,i}=logmean;
            Run.std_mvs{Nlyr,i}=logstd;
            
            
%             Run.md_sig{Nlyr,i}=mean(Run.chains_sig{Nlyr,i}(:,iBurn),2);
%             Run.std_sig{Nlyr,i}=std(Run.chains_sig{Nlyr,i}(:,iBurn),0,2);
            [logmean,logstd]=calc_lognormal_mean(Run.chains_sig{Nlyr,i}(:,iBurn),[Run.GndSigMinLim,Run.GndSigMaxLim]);
            Run.md_sig{Nlyr,i}=logmean;
            Run.std_sig{Nlyr,i}=logstd;
            
            
%             Run.md_sd{Nlyr,i}=mean(Run.chains_sd{Nlyr,i}(:,iBurn),2);
%             Run.std_sd{Nlyr,i}=std(Run.chains_sd{Nlyr,i}(:,iBurn),0,2);
            [logmean,logstd]=calc_lognormal_mean(Run.chains_sd{Nlyr,i}(:,iBurn),[NaN,NaN]);
            Run.md_sd{Nlyr,i}=logmean;
            Run.std_sd{Nlyr,i}=logstd;
            
            
%             Run.md_swe{Nlyr,i}=mean(Run.chains_swe{Nlyr,i}(:,iBurn),2);
%             Run.std_swe{Nlyr,i}=std(Run.chains_swe{Nlyr,i}(:,iBurn),0,2);
            [logmean,logstd]=calc_lognormal_mean(Run.chains_swe{Nlyr,i}(:,iBurn),[NaN,NaN]);
            Run.md_swe{Nlyr,i}=logmean;
            Run.std_swe{Nlyr,i}=logstd;
        end
    end
end



function Run=ModelSelection(Run)   %choose 1,2,3 or N layer plans for each pit?
    Nchan=Run.Np*Run.Nf;
    Cy=Run.StdTb^2*eye(Nchan);
    
    %Note: this hyper-parameter should be read in from file
    switch Run.SturmClass
        case 1
            lambda=2;
            %disp('lambda=5!!!');
        case 2
            lambda=5;
            disp('lambda=5');
        otherwise
            lambda=2;      
    end

    iBurn=Run.Nburn+1:Run.Niter;
    for i=1:Run.Npits,
        Y=Run.TbObs(:,:,i); %TbObs shape is: [Nobs,Nchan,Npit]
        for Nlyr=1:Run.Max_Lyr
            n=Nlyr;
            %Likelihood function. TbPost shape is: [Nchan,Niter,Npits,Nlyr]
            Yhat=median(Run.TbPost(:,iBurn,i,n),2)'; 
            fll(n)=sum(log(mvnpdf(Y,Yhat,Cy)));
            %Number of layers
            p1ll(n)=log(poisspdf(n,lambda));
            %Depth
            dzhat{n}=Run.md_dz{n,i};
            p2ll(n)=log(mvnpdf(log(dzhat{n}),Run.DzMu{n},diag(Run.DzCov{n})));
            %Grain size
            Dhat{n}=Run.md_D{n,i};
            p3ll(n)=log(mvnpdf(log(Dhat{n}),Run.DMu{n},diag(Run.DCov{n})));
            %Density
            Rhohat{n}=Run.md_rho{n,i};
            p4ll(n)=log(mvnpdf(log(Rhohat{n}),Run.RhoMu{n},diag(Run.RhoCov{n})));        
            %Temperature
            That{n}=[Run.md_T{n,i};Run.md_Tsoil{n,i}];
            p5ll(n)=log(mvnpdf(log(That{n}),Run.TMu{n},diag(Run.TCov{n})));                    
            if Run.UsePrior
                Run.J(i,n)=fll(n)+p1ll(n)+p2ll(n)+p3ll(n)+p4ll(n)+p5ll(n);
%                 Run.J(i,n)=fll(n)      +p2ll(n)+p3ll(n)+p4ll(n)+p5ll(n); %the
                %option neglecting the number of layer prior
            else
                Run.J(i,n)=fll(n);
            end
        end
       %disp(num2str(n));   
       %[fll(n),p1ll(n),p2ll(n),p3ll(n),p4ll(n),p5ll(n),p6ll(n),p7ll(n)]'     

    
       switch Run.lyrplan
           case 1
               [~,Run.nHat(i)]=max(Run.J(i,:));
           case 2
               Run.nHat(i)=Run.nlyr_choose;
           case 99
               [~,Run.nHat(i)]=max(Run.J(i,:));
               if(Run.nHat(i)==1)
                   Run.nHat(i)=2;
                   disp('change nhat to 2 from 1:')
                   disp(num2str(Run.nHat(i)))
               end
           otherwise
                disp('please choose how to determin number of layers');
       end
       

        %calculate the final output
        n=Run.nHat(i);
        Run.sdHat(i)=Run.md_sd{n,i};
        Run.SWEHat(i)=Run.md_swe{n,i};
        Run.sdStdHat(i)=Run.std_sd{n,i};
        Run.SWEStdHat(i)=Run.std_swe{n,i};
        Run.TsoilHat(i)=Run.md_Tsoil{n,i};
        Run.MvSHat(i)=Run.md_mvs{n,i};
        Run.SigHat(i)=Run.md_sig{n,i};
        
        mass=Run.md_dz{n,i}.*Run.md_rho{n,i};
        weight=mass/sum(mass);
        
        Run.rhoavgHat(i)=sum(weight.*Run.md_rho{n,i});
        Run.DavgHat(i)=sum(weight.*Run.md_D{n,i});
        Run.TavgHat(i)=sum(weight.*Run.md_T{n,i});
    end
end         



%% plot the result...
function PlotAcceptance(Run,MaxLyr)
    figure
    for i=1:MaxLyr
        subplot(2,MaxLyr/2,i)     
        hist(Run.acceptance.dz(:,i))
        title(['Layer thickness acceptance, Nlyr=' num2str(i) ...
            ', Avg. acpt. =' num2str(mean(Run.acceptance.dz(:,i))) ])           
        hold on;
        a=axis;
        plot([.25 .25],[a(3) a(4)],'r--','LineWidth',2)
        plot([.5 .5],[a(3) a(4)],'r--','LineWidth',2)
        hold off
        set(gca,'xlim',[0,1]);
        xlabel('Acceptance')
        ylabel('Frequency')                
    end
    figure
    for i=1:MaxLyr,
        subplot(2,MaxLyr/2,i)
        hist(Run.acceptance.rho(:,i))
        title(['Density acceptance, Nlyr=' num2str(i) ...
            ', Avg. acpt. =' num2str(mean(Run.acceptance.rho(:,i))) ])                               
        hold on;
        a=axis;
        plot([.25 .25],[a(3) a(4)],'r--','LineWidth',2)
        plot([.5 .5],[a(3) a(4)],'r--','LineWidth',2)
        hold off
        set(gca,'xlim',[0,1]);
        xlabel('Acceptance')
        ylabel('Frequency')                
    end    
    figure
    for i=1:MaxLyr,
        subplot(2,MaxLyr/2,i)
        hist(Run.acceptance.D(:,i))
        title(['Grain size acceptance, Nlyr=' num2str(i) ...
            ', Avg. acpt. =' num2str(mean(Run.acceptance.D(:,i))) ])                                           
        hold on;
        a=axis;
        plot([.25 .25],[a(3) a(4)],'r--','LineWidth',2)
        plot([.5 .5],[a(3) a(4)],'r--','LineWidth',2)
        hold off
        set(gca,'xlim',[0,1]);
        xlabel('Acceptance')
        ylabel('Frequency')                
    end  
    figure
    for i=1:MaxLyr,
        subplot(2,MaxLyr/2,i)
        hist(Run.acceptance.T(:,i))
        title(['Temperature acceptance, Nlyr=' num2str(i) ...
            ', Avg. acpt. =' num2str(mean(Run.acceptance.T(:,i))) ])                                           
        hold on;
        a=axis;
        plot([.25 .25],[a(3) a(4)],'r--','LineWidth',2)
        plot([.5 .5],[a(3) a(4)],'r--','LineWidth',2)
        hold off
        set(gca,'xlim',[0,1]);
        xlabel('Acceptance')
        ylabel('Frequency')                
    end  
end %function PlotAcceptance



function PlotThetaOut(Run,Pit,MaxLyr)
    %X=log([1:Run.Niter]')*2+1;
    %X0=log(Run.Nburn)*2+1;
    X=1:Run.Niter;
    X0=Run.Nburn;
    
    figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        plot(X,Run.chains_dz{j,Pit});
        hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        title(['Pit #' num2str(Pit) ': Layer thickness, Nlyr=' num2str(j)])
        xlabel('Iteration #')
        ylabel('Layer thickness, m')
    end

    figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        plot(X,Run.chains_rho{j,Pit});
        hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        title(['Pit# ' num2str(Pit) ': Layer density, Nlyr=' num2str(j)])
        xlabel('Iteration #')
        ylabel('Density, kg/m^3')                
    end

    figure
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        plot(X,Run.chains_D{j,Pit});
        hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        title(['Pit# ' num2str(Pit) ': Layer grain size, Nlyr=' num2str(j)])
        xlabel('Iteration #')
        ylabel('Grain size, mm')                                
    end
%     
%     figure
%     for j=1:MaxLyr,
%         subplot(MaxLyr/2,2,j)
%         plot(X,274-Run.chains_T{j,Pit}); hold on; 
%         plot(X,274-Run.chains_Tsoil{j,Pit});
%         a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
%         title(['Pit# ' num2str(Pit) ': Layer temperature, Nlyr=' num2str(j)])
%         xlabel('Iteration #')
%         ylabel('Temperature, K')                                
%     end  
%     
%     figure;
%     for j=1:MaxLyr,
%         subplot(MaxLyr/2,2,j)
%         plot(X,Run.chains_mvs{j,Pit});
%         hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
%         title(['Pit# ' num2str(Pit) ': Soil water content, Nlyr=' num2str(j)])
%         xlabel('Iteration #')
%         ylabel('Soil water content, FRAC')                                
%     end  
% 
% 	figure;
%     for j=1:MaxLyr,
%         subplot(MaxLyr/2,2,j)
%         plot(X,Run.chains_sig{j,Pit});
%         hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
%         title(['Pit# ' num2str(Pit) ': Soil roughness, Nlyr=' num2str(j)])
%         xlabel('Iteration #')
%         ylabel('Soil roughness, mm')                                
%     end
%     
%     
%     
% 	figure;
%     for j=1:MaxLyr,
%         subplot(MaxLyr/2,2,j)
%         plot(X,Run.chains_sd{j,Pit});
%         hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
%         title(['Pit# ' num2str(Pit) ': SD, Nlyr=' num2str(j)])
%         xlabel('Iteration #')
%         ylabel('SD, m')                                
%     end
% %     
%     
	figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        plot(X,Run.chains_swe{j,Pit});
        hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        title(['Pit# ' num2str(Pit) ': SWE, Nlyr=' num2str(j)])
        xlabel('Iteration #')
        ylabel('SWE, mm')                                
    end
end %function PlotThetaOut





function PlotThetaOutAvg(Run,Pit,MaxLyr)
    %plot the cumilatied average of the chains
    %X=log([1:Run.Niter]')*2+1;
    %X0=log(Run.Nburn)*2+1;
    X=1:Run.Niter;
    X0=Run.Nburn;
    
    figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        Y=Run.chains_dz{j,Pit};
        Yavg=calc_chains_cumilated_avg(Y,Run.Nburn);
        plot(X,Yavg);
        hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        title(['Pit #' num2str(Pit) ': Layer thickness, Nlyr=' num2str(j)])
        xlabel('Iteration #')
        ylabel('Layer thickness, m')
    end

    figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        Y=Run.chains_rho{j,Pit};
        Yavg=calc_chains_cumilated_avg(Y,Run.Nburn);
        plot(X,Yavg);
        hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        title(['Pit# ' num2str(Pit) ': Layer density, Nlyr=' num2str(j)])
        xlabel('Iteration #')
        ylabel('Density, kg/m^3')                
    end

    figure
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        Y=Run.chains_D{j,Pit};
        Yavg=calc_chains_cumilated_avg(Y,Run.Nburn);
        plot(X,Yavg);
        hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        title(['Pit# ' num2str(Pit) ': Layer grain size, Nlyr=' num2str(j)])
        xlabel('Iteration #')
        ylabel('Grain size, mm')                                
    end
    
    figure
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        
        Y1=274-Run.chains_T{j,Pit};
        Yavg1=calc_chains_cumilated_avg(Y1,Run.Nburn);
        
        Y2=274-Run.chains_Tsoil{j,Pit};
        Yavg2=calc_chains_cumilated_avg(Y2,Run.Nburn);

        plot(X,Yavg1); hold on; 
        plot(X,Yavg2);
        a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        title(['Pit# ' num2str(Pit) ': Layer temperature, Nlyr=' num2str(j)])
        xlabel('Iteration #')
        ylabel('Temperature, K')                                
    end  
    
    figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        Y=Run.chains_mvs{j,Pit};
        Yavg=calc_chains_cumilated_avg(Y,Run.Nburn);
        plot(X,Yavg);
        hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        title(['Pit# ' num2str(Pit) ': Soil water content, Nlyr=' num2str(j)])
        xlabel('Iteration #')
        ylabel('Soil water content, FRAC')                                
    end  

	figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        Y=Run.chains_sig{j,Pit};
        Yavg=calc_chains_cumilated_avg(Y,Run.Nburn);
        plot(X,Yavg);
        hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        title(['Pit# ' num2str(Pit) ': Soil roughness, Nlyr=' num2str(j)])
        xlabel('Iteration #')
        ylabel('Soil roughness, mm')                                
    end
    
	figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        Y=Run.chains_sd{j,Pit};
        Yavg=calc_chains_cumilated_avg(Y,Run.Nburn);
        plot(X,Yavg);
        hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        title(['Pit# ' num2str(Pit) ': SD, Nlyr=' num2str(j)])
        xlabel('Iteration #')
        ylabel('SD, m')                                
    end
    
    
    
	figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        Y=Run.chains_swe{j,Pit};
        Yavg=calc_chains_cumilated_avg(Y,Run.Nburn);
        plot(X,Yavg);
        hold on; a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        title(['Pit# ' num2str(Pit) ': SWE, Nlyr=' num2str(j)])
        xlabel('Iteration #')
        ylabel('SWE, mm')                                
    end
end %function PlotThetaOutAvg






function PlotThetaOutHist(Run,Pit,MaxLyr)
	iBurn=Run.Nburn+1:Run.Niter;
    
    figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(Run.chains_dz{j,Pit}(:,iBurn)',50)
        title(['Pit #' num2str(Pit) ': Layer thickness, Nlyr=' num2str(j)])
        xlabel('Layer thickness, m')
    end

    figure
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(Run.chains_rho{j,Pit}(:,iBurn)',50)
        title(['Pit# ' num2str(Pit) ': Layer density, Nlyr=' num2str(j)])
        xlabel('Density, kg/m^3')                
    end

    figure
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(Run.chains_D{j,Pit}(:,iBurn)',50)
        title(['Pit# ' num2str(Pit) ': Layer grain size, Nlyr=' num2str(j)])
        xlabel('Grain size, mm')                                
    end
    
    figure
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(Run.chains_T{j,Pit}(:,iBurn)',50); hold on;
        hist(Run.chains_Tsoil{j,Pit}(iBurn)',50)
        title(['Pit# ' num2str(Pit) ': Layer temperature, Nlyr=' num2str(j)])
        xlabel('Temperature, K')                                
    end
    
    figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(Run.chains_mvs{j,Pit}(:,iBurn)',50)
        title(['Pit# ' num2str(Pit) ': Soil water content, Nlyr=' num2str(j)])
        xlabel('Soil water content, FRAC')                                
    end  
    
	figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(Run.chains_sig{j,Pit}(:,iBurn)*1000',50)
        title(['Pit# ' num2str(Pit) ': Soil roughness, Nlyr=' num2str(j)])
        xlabel('Soil roughness, mm')                                
    end  
    
    
	figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(Run.chains_sd{j,Pit}(:,iBurn)',50)
        title(['Pit# ' num2str(Pit) ': SD, Nlyr=' num2str(j)])
        xlabel('SD, m')                                
    end  
    
	figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(Run.chains_swe{j,Pit}(:,iBurn)',50)
        title(['Pit# ' num2str(Pit) ': SWE, Nlyr=' num2str(j)])
        xlabel('SWE, mm')                                
    end 
    
    
end %function PlotThetaOut



function PlotThetaLogHist(Run,Pit,MaxLyr)
	iBurn=Run.Nburn+1:Run.Niter;
    
    figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(log(Run.chains_dz{j,Pit}(:,iBurn)'),50)
        title(['Pit #' num2str(Pit) ': Layer thickness, Nlyr=' num2str(j)])
        xlabel('Layer thickness, m')
    end

    figure
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(log(Run.chains_rho{j,Pit}(:,iBurn)'),50)
        title(['Pit# ' num2str(Pit) ': Layer density, Nlyr=' num2str(j)])
        xlabel('Density, kg/m^3')                
    end

    figure
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(log(Run.chains_D{j,Pit}(:,iBurn)'),50)
        title(['Pit# ' num2str(Pit) ': Layer grain size, Nlyr=' num2str(j)])
        xlabel('Grain size, mm')                                
    end
    
    figure
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        %hist(274-Run.chains_T{j,Pit}(:,iBurn)',50)
        %hist(274-Run.chains_Tsoil{j,Pit}(iBurn)',50)
        hist(log(Run.chains_T{j,Pit}(:,iBurn)'),50);hold on;
        hist(log(Run.chains_Tsoil{j,Pit}(iBurn)'),50)
        title(['Pit# ' num2str(Pit) ': Layer temperature, Nlyr=' num2str(j)])
        xlabel('Temperature, K')                                
    end
    
    figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(log(Run.chains_mvs{j,Pit}(:,iBurn)'),50)
        title(['Pit# ' num2str(Pit) ': Soil water content, Nlyr=' num2str(j)])
        xlabel('Soil water content, FRAC')                                
    end  
    
	figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(log(Run.chains_sig{j,Pit}(:,iBurn)*1000'),50)
        title(['Pit# ' num2str(Pit) ': Soil roughness, Nlyr=' num2str(j)])
        xlabel('Soil roughness, mm')                                
    end  
    
    
	figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(log(Run.chains_sd{j,Pit}(:,iBurn)'),50)
        title(['Pit# ' num2str(Pit) ': SD, Nlyr=' num2str(j)])
        xlabel('SD, m')                                
    end  
    
	figure;
    for j=1:MaxLyr,
        subplot(MaxLyr/2,2,j)
        hist(log(Run.chains_swe{j,Pit}(:,iBurn)'),50)
        title(['Pit# ' num2str(Pit) ': SWE, Nlyr=' num2str(j)])
        xlabel('SWE, mm')                                
    end 
    
    
end %function PlotThetaLogHist





function PlotThetaCompare(Run,Pit)

    %the MCMC result
    Nlyr=Run.nHat(Pit);
    
    dz_post=Run.md_dz{Nlyr,Pit};    dz_post=squeeze(dz_post); %this is estimated snow thickness
    %
%     dz_post=dz_post*0.76; %revised for Lidar test only
    
    rho_post=Run.md_rho{Nlyr,Pit};  dz_post=squeeze(dz_post); 
    D_post=Run.md_D{Nlyr,Pit};      dz_post=squeeze(dz_post); 
    T_post=Run.md_T{Nlyr,Pit};      dz_post=squeeze(dz_post); 
    
    %the truth
%     Pit=8; %revised for Lidar test only
    dz_true=Run.TruePits(Pit).dz;%to be revised
    rho_true=Run.TruePits(Pit).density; %to be revised
    T_true=Run.TruePits(Pit).T;  %to be revised
    if(Run.Obs_Model==1)
        D_true=Run.TruePits(Pit).pex;  %to be revised
    else
        D_true=Run.TruePits(Pit).dmax;  %to be revised
    end

    %plot and compare
    figure;
    set(gcf,'color','w');
    subplot(1,3,1);
    plot_thisway(Run,dz_post*100,rho_post,'ro-');
    plot_thisway(Run,dz_true*100,rho_true,'ko-');
    title(['pits', num2str(Pit), ': density profile']);
    xlabel('density (kg/m^3)'); ylabel('z (cm)'); legend('post','true');
    
    subplot(1,3,2);
    plot_thisway(Run,dz_post*100,D_post,'ro-');
    plot_thisway(Run,dz_true*100,D_true,'ko-');
    title('Dmax profile');
    xlabel('Dmax (mm)'); ylabel('z (cm)');legend('post','true');
    
    subplot(1,3,3);
    plot_thisway(Run,dz_post*100,274-T_post-273.15,'ro-');
    plot_thisway(Run,dz_true*100,T_true-273.15,'ko-');
    title('Temp profile')
    xlabel('T (^oC)'); ylabel('z (cm)');legend('post','true');                        
end



function plot_thisway(Run,z,data,Mark)
    depth=sum(z);
    z2=0.1:0.1:depth;
    
    data2=zeros(0,1);
    for k=1:length(z)
        data2=[data2;repmat(data(k),round(z(k)*10),1)];
    end
    
    z2=z2-depth;
    
    if(length(z2)>length(data2))
        disp(['length(z2)>length(data2) by ', num2str(length(z2)-length(data2))])
        z2=z2(1:length(data2));
    elseif (length(z2)<length(data2))
        disp(['length(z2)<length(data2) by ', num2str(length(z2)-length(data2))])
        data2=data2(1:length(z2));
    end
    plot(data2,z2,Mark); hold on;
end



end %methods

end %classdef
