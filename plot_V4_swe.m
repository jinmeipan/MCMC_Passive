%plot_V4_results_general
%Nov28, 2015, add percentile calculation

model='iba';


Strings={'sd_nov','sd_dec','sd_jan','sd_feb','sd_mar','ch_jan','ch_feb','ch_mar','ch_apr','clpx_feb'};
labelstr={'mo','m+','m^','msq','md','bo','b*','b^','bsq','kd'};
planStr={'gSWE-bStra-real','bSWE-bStra-real','gSWE-gStra-real','bSWE-gStra-real',...
    'gSWE-bStra-syns','bSWE-bStra-syns','gSWE-gStra-syns','bSWE-gStra-syns'};

%
planc=1:8;
sitemonc=1:5; %10;


% figure;
set(gcf,'color','w');

folder='/Users/office/Downloads/MCMC/MCMCRunData_V2/';
fidn=fopen([folder,'!!swe_stat_',model,'.txt'],'w');
fprintf(['SWE RMSE estimate','SWE RRMSE estimate','SWE bias estimate',...
    'SWE RMSE prior','SWE RRMSE prior','SWE bias prior']);


% iz=[2,1,4,3]+4;
iz=2;
% iz=1:1:8;

for ik=1:length(iz) %1:length(planc) 
    
    i=iz(ik);
    subplot(4,2,ik)
    
    SWE_MCMC=nan;
    SWE_TRUE=nan;
    Nnn=nan;
    SWE_Percentile=nan;
    SWE_Prob=nan;
	Np=nan(length(planc),1);
    
    
    for j=1:length(sitemonc)
        
        rfile_option='rfile3-test2freq';
        
        
        name=['mcmc',num2str(i),'-',model,num2str(j)]
        load([folder,rfile_option,'/',name,'.mat'],'r');
        swe_mcmc = r.SWEHat;
        
        
%         exclude two pits using IBA
%         if(strcmp(model,'iba')==1)
%             if(j==8 | j==9) % | j==10)
%                 disp(swe_mcmc);
%                 swe_mcmc=swe_mcmc*NaN;
%             end
%         end
        
%         if(strcmp(model,'iba')==1)
%             if(j==8 | j==9 | j==10)
%                 disp(swe_mcmc);
%                 swe_mcmc=swe_mcmc*NaN;
%             end
%         end
%         
                
        
        disp(['plan-sitem:',num2str([i,j])])
        disp(num2str(r.nHat))
        
        swe_true=nan(r.Npits,1);
        Nnn=[Nnn;r.nHat'];
        perct=nan(r.Npits,1);
        prob=nan(r.Npits,1);
        
        for k=1:r.Npits
            swe_true(k)=r.TruePits(k).SWE;
            
            %percentile
            swe_chains=r.chains_swe{r.nHat(k),k};
            swe_chains=swe_chains(r.Nburn:r.Niter);
            index=find(swe_chains<=0);
            if(length(index)>0)
                disp('swe less than zero!')
                disp(num2str(swe_chains(index)))
                swe_chains(index)=[];
            end
            
            %percentile of trueswe
            if(0)
                swe_t=swe_true(k);
                mean_log=mean(log(swe_chains));
                std_log= std(log(swe_chains)); %sqrt(sum((log(swe_chains)-mean_log).^2)/ (length(swe_chains)-1));

                log_swe_t=log(swe_t);       
                perct(k)=(log_swe_t-mean_log)./std_log;
                prob(k)= 1/(std_log*sqrt(2*pi())) * exp (( - (log_swe_t-mean_log).^2 )/ (2*std_log.^2));
            else
                prob(k)=nan;
                
                swe_t=swe_true(k);
                swe_r=r.md_swe{r.nHat(k),k};
                
                index=find(swe_chains==swe_t);
                if(length(index)>0)
                    error('same values for swe_t!');
                end
                
                index=find(swe_chains==swe_r);
                if(length(index)>0)
                    error('same values for swe_r!');
                end
                
                [B,IX]=sort([swe_t,swe_chains]);
                [B2,IX2]=sort([swe_r,swe_chains]);
                
                index1=find(B==swe_t);
                index2=find(B2==swe_r);
                
                perct(k)=(index1-index2)/(length(swe_chains)+1);
            end
            
            
            
            clear swe_t swe_r check_samevalue B IX swe_chains
            clear index1 index2 index B2 IX2
        end
        
        
        %better!!
%         if(i==2 & j==5)
%             clear r          
%             load([folder,rfile_option,'/mcmc2-iba5-pit4.mat'],'r');
%             swe_mcmc(4)=r.SWEHat(1);
%         end
%         
%         if(i==8 & j==2) 
%             clear r
%             load([folder,rfile_option,'/mcmc8-iba2-pit1.mat'],'r');
%             swe_mcmc(1)=r.SWEHat(1);
%         end
        
        
        plot(swe_true,swe_mcmc,labelstr{j}); hold on;
        clear r
        
        SWE_MCMC=[SWE_MCMC;swe_mcmc'];
        SWE_TRUE=[SWE_TRUE;swe_true];
        SWE_Percentile=[SWE_Percentile;perct];
        SWE_Prob=[SWE_Prob;prob];
        Np(j)=length(swe_true);
    end
    
    if(i==1)
        Legendstr={'sd-nov','sd-dec','sd-jan','sd-feb','sd-mar','ch-jan','ch-feb','ch-mar','ch-apr','clpx-feb'};
        legend(Legendstr,'location','northwest');
    end
    
    
    plot([0,400],[0,400],'k--');
    title(planStr{i});
    xlabel('SWE_{true} (mm)');
    ylabel('SWE_{MCMC}(mm)')
    
    
    %calculate error!
    SWE_MCMC(1)=[];
    SWE_TRUE(1)=[];
    Nnn(1)=[];
    SWE_Percentile(1)=[];
    SWE_Prob(1)=[];
	[SWE_PR,SiteNo]=get_pr_swe(planStr{i},Strings,Np);
    
    index=find(isnan(SWE_MCMC)~=1);
    SWE_MCMC=SWE_MCMC(index);
    SWE_TRUE=SWE_TRUE(index);
    Nnn=Nnn(index);
    SWE_PR=SWE_PR(index);
    SiteNo=SiteNo(index);
    SWE_Percentile=SWE_Percentile(index);
    SWE_Prob=SWE_Prob(index);
    
    
    %mcmc error
	N=length(SWE_MCMC);
    rms_mcmc= sqrt(sum((SWE_MCMC-SWE_TRUE).^2)/(N-1));
    rrms_mcmc = sqrt(sum( ((SWE_MCMC-SWE_TRUE)./SWE_TRUE).^2)/(N-1));
    bias_mcmc = mean(SWE_MCMC-SWE_TRUE);
    perct_mcmc = mean(abs(SWE_Percentile));
    prob_mcmc = mean(SWE_Prob);
    
    %prior error
    rms_pr= sqrt(sum((SWE_PR-SWE_TRUE).^2)/(N-1));
	rrms_pr = sqrt(sum( ((SWE_PR-SWE_TRUE)./SWE_TRUE).^2)/(N-1));
    bias_pr = mean(SWE_PR-SWE_TRUE);
    
    
    fprintf(fidn,[num2str([i, rms_mcmc, rrms_mcmc, bias_mcmc, rms_pr, rrms_pr, bias_pr, perct_mcmc, prob_mcmc]), '\n']);
    
    

    fidn2=fopen([folder,'!!swe_stat_',model,'_plan',num2str(i),'.txt'],'w');
    fprintf(fidn2, ['SWE_TRUE   SWE_PR    SWE_MCMC  N SWE_Percentile  SWE_prob \n']);
    
    for ik=1:length(SWE_TRUE)
        fprintf(fidn2,[num2str([SiteNo(ik),SWE_TRUE(ik),SWE_PR(ik),SWE_MCMC(ik),Nnn(ik),SWE_Percentile(ik),SWE_Prob(ik)]),'\n']);
    end
    fclose(fidn2);
  
end

fclose(fidn);


Mtx=[SWE_TRUE,SWE_MCMC,Nnn];
open Mtx
