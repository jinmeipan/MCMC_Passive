%inputs:
%-sp: in the format of snowpit structure; e.g. sp=tsp_sdkl_nov;
%-site: where is the snowpits, to determine the frequencies.
%-stdTb: observation error of TB (K); eg. stdTb=2.0



function prep_tbobs(sp,freq,freq_index,theta,tb_type,stdTb,site,model,folder0,file,TskyStr)

disp(['number of pits=',num2str(length(sp))])
disp('check observed tb frequencies to be chosen=')


%import true tb anyway
for i=1:length(sp)
    switch site
        case 'sdkl'
            TB_true(i,:)=sp(i).tbv(freq_index);
        case 'churchill'
            TB_true(i,:)=sp(i).tbv(freq_index);
        case 'clpx'
            TB_true(i,:)=sp(i).tbv(freq_index);
        otherwise
            error('please determine which frequency to extract');
    end
end


switch tb_type
    case 'real'
        TB_MCMC=TB_true;
        calc_syns_tb=0;
    case 'syns'
        %synthetic tb
        calc_syns_tb=1;
    otherwise
        error('please choose a type of tb to be inputted to mcmc: real or syns!');
end



%%
if(calc_syns_tb==1)
    %synthetic tb
	TB_MCMC0=nan(length(sp),length(freq));
    for i=1:length(sp);
        %frequencies, theta
        switch site
            case 'sdkl'
                gndsig=0.3/100; %cm to m; %better change? from 0 to 0.3, to prevent non-zero estimate
                typ_eps_fix=0;
                eps_soil=0;
                mv_soil=2; %percent

            case 'churchill'       
                gndsig=0.4/100; %cm to m;
                typ_eps_fix=0;
                eps_soil=0;
                mv_soil=12; %percent
 
            case 'clpx'
                gndsig=0.1/100; %cm to m;
                typ_eps_fix=0;
                eps_soil=0;
                mv_soil=12; %percent
        end
        
        
        %model-specific
        rho=1.5;
        vsand= 30.6299992;
        vsilt=55.8899994;
        vclay= 13.4799995;
        %Tsky=2.7; %jinmei, revised Dec22,2015
        
        switch model
            case 'iba'
                layernumber=[1:sp(i).nlayer]';
                thickness=sp(i).dz*100;  %cm
                density=sp(i).density;
                T=sp(i).T;
                Zeros=zeros(sp(i).nlayer,1);
                pex=sp(i).pex;
                Tgnd=sp(i).soilT;
                Tsoil=Tgnd-273.15;
%                 mv_soil=sp(i).mv_soil*100;
                scopt=1; %1=ise born approximation
                
                
                y=[layernumber,T,Zeros,density,thickness,Zeros,pex];
                for j=1:length(freq)
                    
                    Tsky=TskyStr(j); %jinmei added, Dec22,2015
                            
                    f=freq(j);
                    [s0v,s0h]=soil_reflectivity2(f,theta,y,Tsoil,mv_soil,gndsig,rho,vsand,vsilt,vclay,...
                            typ_eps_fix,eps_soil);
                    [Tbh,Tbv] = memlsfunc(f,theta,s0h,s0v,y,Tsky,Tgnd,scopt);
                    TB_MCMC0(i,j)=Tbv;
                end
                clear layernumber thickness density T Zeros pex Tgnd Tsoil scopt y f j s0v s0h Tbh Tbv
                
            case 'emp'
                layernumber=[1:sp(i).nlayer]';
                thickness=sp(i).dz*100;  %cm
                density=sp(i).density;
                T=sp(i).T;
                Zeros=zeros(sp(i).nlayer,1);
                pex=sp(i).pex;
                Tgnd=sp(i).soilT;
                Tsoil=Tgnd-273.15;
%                 mv_soil=sp(i).mv_soil*100;
                scopt=2; %2=empirical memls
                
                
                y=[layernumber,T,Zeros,density,thickness,Zeros,pex];
                for j=1:length(freq)
                   
                    Tsky=TskyStr(j); %jinmei added, Dec22,2015
                    
                    f=freq(j);
                    [s0v,s0h]=soil_reflectivity2(f,theta,y,Tsoil,mv_soil,gndsig,rho,vsand,vsilt,vclay,...
                            typ_eps_fix,eps_soil);
                    [Tbh,Tbv] = memlsfunc(f,theta,s0h,s0v,y,Tsky,Tgnd,scopt);
                    TB_MCMC0(i,j)=Tbv;
                end
                clear layernumber thickness density T Zeros pex Tgnd Tsoil scopt y f j s0v s0h Tbh Tbv
                                
                                
            case 'hut' %here we are thinking to use the combined model. 
                layernumber=sp(i).nlayer;
                thickness=sp(i).dz;  %m
                density=sp(i).density;
                T=sp(i).T-273.15;
                Zeros=zeros(layernumber,1);
                dmax=sp(i).dmax;
                deff=dmax2deff(dmax);
                Tsoil=sp(i).soilT-273.15;
%                 mv_soil=sp(i).mv_soil*100;
                
                for j=1:length(freq)
                    
                    Tsky=TskyStr(j); %jinmei added, Dec22,2015
                    
                    f=freq(j); 
                    if(f<30.0)
                        typke=1; %1=Hallikainen eq.
                    else
                        typke=2; %2=Roy eq.
                    end
                    
                    [TB_v,TB_h,rss_v,rss_h]=multi_HUT(f, theta,...
                            layernumber, thickness,T,density,Zeros, deff,...
                            gndsig,Tsoil,mv_soil,rho,vsand,vsilt,vclay,...
                            Tsky,Tsky,...
                            typke,typ_eps_fix);
                    TB_MCMC0(i,j)=TB_v;
                end
                clear layernumber thickness dentisy T Zeros dmax deff Tsoil j f typke 
                clear TB_v TB_h rss_v rss_h
        end
    end
    
    %add random error
    TB_MCMC=TB_MCMC0+stdTb*randn(length(sp),length(freq));
    
% 	for i=1:length(sp)
%         disp([file,':'])
%         disp(['pit ',num2str(i)])
%         disp(num2str(TB_true(i,:)));
%         disp(num2str(TB_MCMC0(i,:)));
%         disp(num2str(TB_MCMC(i,:)));
%         disp(' ')
%     end

end


%% write
fid_tb=fopen([folder0,'temp/tb_obs_',file,'.txt'],'w');
for ip=1:size(TB_MCMC,1)   
    fprintf(fid_tb,[' Pit  #     '     num2str(ip) '  \n']);

    for j=1:size(TB_MCMC,2)
        fprintf(fid_tb,'%6.1f\n', TB_MCMC(ip,j));
    end
end
fclose(fid_tb);

clear TB_MCMC TB_true TB_MCMC0
clear ip j i fid_tb
clear gndsig typ_eps_fix eps_soil mv_soil rho vsand vsilt vclay Tsky

end
