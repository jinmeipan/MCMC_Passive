% All MCMC fortran code input files preparation
% this is for all R3 results
clear all

folder0='/Users/office/Downloads/MCMC/MCMCRunData_V2/MCMC_4FREQ_NEW/MCMC_2/';
mkdir(folder0);
cd(folder0)
if(0)
    open prep_filename
    open prep_runparam2
    open prep_tbobs
    open prep_hyper_gp
    open prep_hyper_gp_mlyr
    open prep_hyper_bp
    open prep_hyper_bp_mlyr
    open prep_truesp
end


%% Basic settings !!!
prior_type='bp';  %'gp'=good prior; 'bp'=vic or nldas prior
stra_type='nostra';  %'stra'=use strategraphy; 'i=nostra'=use homogenous prior
tb_type='real'; %'real'=real observation; 'syns'=synthetic tb with stdTb;


model='iba';   %'emp'=empirical memls; 'iba'='iba','hut'=hut model combination
stdTb=2.0;


freq_plan='4freq'; %'4freq' or '3freq'.

%%
testyes=0;
prep_allfiles(folder0,prior_type,stra_type,tb_type,model,stdTb,testyes,freq_plan);

