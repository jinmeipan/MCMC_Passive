%use Note to exclude the point that hit the max,min_limit

function [logmean,logstd]=calc_lognormal_mean(yy,Note)

%%
% iBurn=r.Nburn+1:r.Niter; 
% yy=r.chains_swe{2,1}(:,iBurn); Note=[0,10000];
% % chains_rho
% % chains_D
% % chains_T
% % chains_Tsoil
% % chains_mvs
% % chains_sig
% % chains_sd
% % chains_swe
% figure;hist(yy'); hold on;
% clc

%type
% type=1 ;%use maximum possibility
% type=2 ;%use direct mean
% type=3 ;%use logmean
% type=4 ;%average of 95% percentile
type=5 ;%use logmean or direct mean, judged by skewness

logmean=nan(size(yy,1),1);
logstd=nan(size(yy,1),1);


for ir=1:size(yy,1)
    
    yy2=yy(ir,:);
    
    %remove values close to the edges
    if(isnan(Note(1))==0 & isnan(Note(2))==0)
        index=find(yy2<=Note(1)*(1+1e-3) | yy2>=Note(2)*(1-1e-3));
        yy2(index)=[];
    end
    
    %
    switch type
        case 1
            %pick the nbin with highest probability
            [n,bin]=hist(yy2,50);
            
            index=find(n==max(n));
            if(length(index)>1)
                disp('there are two nbins have the same largest number of frequency');
                disp(num2str(index))
                disp(num2str(n(index)))
            end
            

            %pick the values fallen into this bin
            if(sum(index>1)==length(index))
                x_low= (bin(index)+bin(index-1))/2;
            else
                x_low=min(yy2);
            end

            if(sum(index<length(bin))==length(index))
                x_high=(bin(index)+bin(index+1))/2;
            else
                x_high=max(yy2);
            end

            if(length(index)>1)
                x_low=min(x_low);
                x_high=max(x_high);
            end

            x_index=find(yy2>x_low & yy2<x_high);

            %take the average of the values fallen into this bin
            %because there are 28000 (or 8000) iterations, there are 560 (or 160) values 
            %in this bin
            logmean(ir)=mean(yy2(x_index));
            
        case 2
            logmean(ir)=mean(yy2);
            
        case 3
            index=find(yy2>0);
            log_yy2=log(yy2(index));
            logmean(ir)=exp(mean(log_yy2));
            
        case 4
            n=length(yy2);
            [yy3,index]=sort(yy2);
            index2=fix(n/2-n*0.025): 1: fix(n/2+n*0.025);
            logmean(ir)=mean(yy3(index2));
        case 5
            index=find(yy2>0);
            log_yy2=log(yy2(index));
            
            sk_log_yy2=skewness(log_yy2);
            sk_yy2=skewness(yy2);
            
            %disp('sk_yy2 sk_log_yy2');
            %disp(num2str([sk_yy2,sk_log_yy2]))

            
            if(abs(sk_log_yy2)<abs(sk_yy2) & sk_yy2>0) %sk_yy2>0 means, direct mean leans to left
                logmean(ir)=exp(mean(log_yy2));
            else
                logmean(ir)=mean(yy2);
            end
    end
    
	%s.t.d.
	dum=(yy2-logmean(ir)).^2;
	logstd(ir)= sqrt(sum(dum)/length(dum-1));
end

%%
end
