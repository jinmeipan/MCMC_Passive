%script to plot the output from the snow MCMC

cd('/Users/office/Downloads/MCMC/MCMCRunData_V2');
%RunExp



%%
cd('/Users/office/Downloads/MCMC/MCMCRunData_V2');
r=MCMCRun2;
r.filename='/Users/office/Downloads/MCMC/MCMCRunData_V2/MCMC_Test2freq/MCMC_2/iba3/FILENAME.txt';
% r.filename='/Users/office/Downloads/MCMC/MCMCRunData_V2/MCMC_IBA/iba5/FILENAME.txt';
% r.filename='/Users/office/Downloads/MCMC/MCMCRunData_V2/MCMC_IBA/4FREQ_CLPX/MCMC_8/iba2/FILENAME.txt';
% r.filename='/Users/office/Downloads/MCMC/MCMCRunData_V2/MCMC_IBA/MCMC_IBA_LC1/MCMC_2/iba2/FILENAME.txt';
% r.filename='/Users/office/Downloads/MCMC/MCMCRunData_V2/MCMC_4FREQ_NEW/MCMC_2/iba9/FILENAME.txt';
r.lyrplan=99; %1=use modelSelection to select number of layers; 99=contrained minimum number of layers to 2; 2=use r.nlyr_choose
r.nlyr_choose=4; %write when lyplan=2
r.run_location=0; %=1: on supper computer,=0 on desktop
r.Max_Lyr=2;
r=r.Main;


% for i=1:r.Npits   
% 	MaxLyr=2; %up to which number of layer plan
% % 	r.PlotThetaOut(i,MaxLyr);
% 	%r.PlotThetaOutAvg(i,MaxLyr);
% 	%r.PlotThetaOutHist(i,MaxLyr);
%     %r.PlotThetaLogHist(i,MaxLyr);
%     r.PlotThetaCompare(i);
% end


% figure;
% for i=1:10
% plot(exp(r.J(8,:)));hold on;
% end

% r.PlotThetaCompare(2);
% 
% r.SWEHat(8)
% r.PlotThetaCompare(8);
%%
if(true)
swe=r.SWEHat(8)';
disp(swe);

%check
Run=r;
Pit=8;
iplot=0; %0 or 3

%the MCMC result
    Nlyr=Run.nHat(Pit);
    
    dz_post=Run.md_dz{Nlyr,Pit};    dz_post=squeeze(dz_post); %this is estimated snow thickness
    rho_post=Run.md_rho{Nlyr,Pit};  dz_post=squeeze(dz_post); 
    D_post=Run.md_D{Nlyr,Pit};      dz_post=squeeze(dz_post); 
    T_post=Run.md_T{Nlyr,Pit};      dz_post=squeeze(dz_post); 
    
    %the truth
    dz_true=Run.TruePits(Pit).dz;%to be revised
    rho_true=Run.TruePits(Pit).density; %to be revised
    T_true=Run.TruePits(Pit).T;  %to be revised
    if(Run.Obs_Model==1)
        D_true=Run.TruePits(Pit).pex;  %to be revised
    else
        D_true=Run.TruePits(Pit).dmax;  %to be revised
    end

    subplot(2,3,iplot+1)
    plot_thisway(Run,dz_post*100,rho_post,'ro-');
    plot_thisway(Run,dz_true*100,rho_true,'ko-');
    title(['pits', num2str(Pit), ': density profile']);
    xlabel('density (kg/m^3)'); ylabel('z (cm)'); legend('post','true');
    

    subplot(2,3,iplot+2)
    plot_thisway(Run,dz_post*100,D_post,'ro-');
    plot_thisway(Run,dz_true*100,D_true,'ko-');
    title('Dmax profile');
    xlabel('Dmax (mm)'); ylabel('z (cm)');legend('post','true');
    

    subplot(2,3,iplot+3)
    plot_thisway(Run,dz_post*100,274-T_post-273.15,'ro-');
    plot_thisway(Run,dz_true*100,T_true-273.15,'ko-');
    title('Temp profile')
    xlabel('T (^oC)'); ylabel('z (cm)');legend('post','true');          
end  
    
%%

% r.PlotAcceptance(r.Max_Lyr);

nlyr=r.nHat';
sd=r.sdHat';
swe=r.SWEHat';
density=r.rhoavgHat';
gs=r.DavgHat';
T=r.TavgHat';
tsoil=r.TsoilHat';
mvsoil=r.MvSHat';
sigsoil=r.SigHat';

M=[sd;nan;nan;swe;nan;nan;density;nan;nan;gs;nan;nan;T;nan;nan;tsoil;nan;nan;...
    mvsoil;nan;nan;sigsoil];
M(1:r.Npits,2)=nlyr;

open M



%save the results in the format of Structure MCMC automatically
%with names to be revised!
string=r.filename_tbobs;
temp=length('/Users/office/Downloads/MCMC/MCMCRunData_V2/!tb_obs_');
string2=string(temp+1:end-4);
temp2=length(string2);
if(temp2>7)
    complexname=1;
    string2=string2(1:end-3);
end

% name=['bp_hut_',string2]   %!!to be modified!

if(r.lyrplan==2)
    name=[name,'_2lyr!'];
end
save(['/Users/office/Downloads/MCMC/MCMCRunData_V2/Z_Results/',name],'r')


%% check theta post
%theta_post: theta_post{Nlyr}=nan(Ntheta,Run.Niter,Run.Npits)
for i=1:r.Npits   
	MaxLyr=2; %up to which number of layer plan
% 	r.PlotThetaOut(i,MaxLyr);
	%r.PlotThetaOutAvg(i,MaxLyr);
% 	r.PlotThetaOutHist(i,MaxLyr);
    %r.PlotThetaLogHist(i,MaxLyr);
    r.PlotThetaCompare(i);
end

%%
r.PlotAcceptance(r.Max_Lyr);

%%
%check tbpost
%TbPost: [Nchan,Niter,Npits,Nlyr]
for npit=1:r.Npits

% nplan=2;
nplan=r.nHat(npit);

no_freq=3;

switch no_freq
    case 3

        figure;
        set(gcf,'Position',[300,300,900,800]);
        plot(r.TbPost(1:3,:,npit,nplan)');hold on;  
        plot(r.Niter,r.TbObs(:,1,npit),'Marker','o','MarkerSize',10,'color','blue')
        plot(r.Niter,r.TbObs(:,2,npit),'Marker','o','MarkerSize',10,'color','green')
        plot(r.Niter,r.TbObs(:,3,npit),'Marker','o','MarkerSize',10,'color','red')
        legend('10.65','36.5','89','10.65O','36.5O','89O','location','north');hold on;
        title(['Pit ',num2str(npit),' ,V pol.'])
    case 4
        figure;
        set(gcf,'Position',[300,300,900,800]);
        plot(r.TbPost(1:4,:,npit,nplan)');hold on;  
        plot(r.Niter,r.TbObs(:,1,npit),'Marker','o','MarkerSize',10,'color','blue')
        plot(r.Niter,r.TbObs(:,2,npit),'Marker','o','MarkerSize',10,'color','green')
        plot(r.Niter,r.TbObs(:,3,npit),'Marker','o','MarkerSize',10,'color','red')
        plot(r.Niter,r.TbObs(:,4,npit),'Marker','o','MarkerSize',10,'color','yellow')
        legend('10.65','18.7','36.5','89','10.65O','18.7O','36.5O','89O','location','north');hold on;
        title(['Pit ',num2str(npit),' ,V pol.'])
    case 2
        figure;
        set(gcf,'Position',[300,300,900,800]);
        plot(r.TbPost(1:2,:,npit,nplan)');hold on;  
        plot(r.Niter,r.TbObs(:,1,npit),'Marker','o','MarkerSize',10,'color','blue')
        plot(r.Niter,r.TbObs(:,2,npit),'Marker','o','MarkerSize',10,'color','green')
        legend('18.7','36.5','location','north');hold on;
        title(['Pit ',num2str(npit),' ,V pol.'])
end




% figure;plot(r.TbPost(5:8,:,npit,nplan)');hold on;  
% plot(10000,r.TbObs(:,5,npit),'Marker','o','MarkerSize',10,'color','blue')
% plot(10000,r.TbObs(:,6,npit),'Marker','o','MarkerSize',10,'color','green')
% plot(10000,r.TbObs(:,7,npit),'Marker','o','MarkerSize',10,'color','red')
% plot(10000,r.TbObs(:,8,npit),'Marker','o','MarkerSize',10,'color','k')
% legend('10.65','18.7','36.5','89','10.65O','18.7O','36.5O','89O','location','north');hold on;
% title('H pol.')
end

%%
npit=8;
figure;
hist(r.TbPost(1:4,iBurn,npit,nplan)',100); hold on;
for i=1:4
    plot([r.TbObs(:,i,npit),r.TbObs(:,i,npit)],[0,2000],'k--');
end


for i=1:4
    plot(r.TbObs(:,i,npit),1800,'kx');
end

%%
%check chains
pit=1;
lyr=2;
iBurn=r.Nburn+1:r.Niter;
figure;plot(r.chains_dz{lyr,pit}(2,iBurn),r.chains_rho{lyr,pit}(2,iBurn),'bo');xlabel('thickness(m)');ylabel('density(kg/m^3)')
figure;plot(r.chains_rho{lyr,pit}(2,iBurn),r.chains_D{lyr,pit}(2,iBurn),'bo');xlabel('density(kg/m^3)');ylabel('pex(mm)')
figure;plot(r.chains_D{lyr,pit}(2,iBurn),r.chains_T{lyr,pit}(2,iBurn),'bo');ylabel('T(degC)');xlabel('pex(mm)')

figure;plot(r.chains_rho{lyr,pit}(1,iBurn),r.chains_rho{lyr,pit}(2,iBurn),'bo');
figure;plot(r.chains_T{lyr,pit}(1,iBurn),r.chains_T{lyr,pit}(2,iBurn),'bo');
