% All MCMC fortran code input files preparation
% this is for all R4 results





function prep_allfiles(folder0,prior_type,stra_type,tb_type,model,stdTb,testyes,freq_plan)


    if(testyes==1)
        site_string='sdkl' %'sdkl'=sodankyla; 'churchill'; 'clpx'
        month_string=3;    %in numbers
    else
        site_string={'sdkl','sdkl','sdkl','sdkl','sdkl','churchill','churchill','churchill','churchill','clpx'};  %????????
        month_string=[11,12,1,2,3,    1,2,3,4,    2];  
    end


%     for ikk=1:10
    for ikk=10:10
        site=site_string{ikk};      %'sdkl'=sodankyla; 'churchill'; 'clpx'
        month=month_string(ikk)     %2;

        %%
        %file
        mkdir([folder0,'/temp/']);

        % sd,ch,clpx+mon
        switch site
            case 'sdkl'
                sitestr='sd';
            case 'churchill'
                sitestr='ch';
            case 'clpx'
                sitestr='clpx';
        end
        switch month
            case 11
                monthstr='nov';
            case 12
                monthstr='dec';
            case 1
                monthstr='jan';
            case 2
                monthstr='feb';
            case 3
                monthstr='mar';
            case 4
                monthstr='apr';
        end
        file=[sitestr,monthstr,model];
        foldernames=[model,num2str(ikk)];


        switch site
            case 'clpx'
                clpx_use_nldas=1; %!!!! to be determined, there is a third choice for clpx bad priors
            otherwise
                clpx_use_nldas=0; %!!!! to be determined, there is a third choice for clpx bad priors
        end


        %% true_snowpits read for prep_tbs
        switch site
            case 'sdkl'
                load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_sdkl.mat');
                switch month
                    case 11
                        sp=tsp_sdkl_nov;
                    case 12
                        sp=tsp_sdkl_dec;
                    case 1
                        sp=tsp_sdkl_jan;
                    case 2
                        sp=tsp_sdkl_feb;
                    case 3
                        sp=tsp_sdkl_mar;
                end
            case 'churchill'
                load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_churchill.mat');
                switch month
                    case 1
                        sp=tsp_churchill_jan;
                    case 2
                        sp=tsp_churchill_feb;
                    case 3
                        sp=tsp_churchill_mar;
                    case 4
                        sp=tsp_churchill_apr;
                end
            case 'clpx'
                load('/Users/office/Desktop/CLPX_data_FromNSIDC/tsp_clpx_lsos.mat');
                sp=sp_clpx_lsos;
        end


        %%
        %run files
        npol=1;

        prep_filename(folder0,file);
        
        [freq_index,freq,theta,TskyStr]=prep_runparam2(sp,model,site,npol,stdTb,folder0,file,freq_plan); %%
        prep_tbobs(sp,freq,freq_index,theta,tb_type,stdTb,site,model,folder0,file,TskyStr); %%
        prep_truesp(sp,model,folder0);
        switch prior_type
            case 'gp'
                switch stra_type
                    case 'stra'
                        prep_hyper_gp_mlyr(site,model,month,folder0,file);
                    case 'nostra'
                        prep_hyper_gp(site,model,month,folder0,file);
                end
            case 'bp'
                switch stra_type
                    case 'stra'
                        prep_hyper_bp_mlyr(site,model,month,folder0,file,clpx_use_nldas);
                    case 'nostra'
                        prep_hyper_bp(site,model,month,folder0,file,clpx_use_nldas);
                end
        end


        %%
        %change folder names
    %     copy_mcmc_program(model,folder0)

    
        copyfile([folder0,'temp/'],[folder0,foldernames])
        rmdir([folder0,'temp/'],'s')

    end

end