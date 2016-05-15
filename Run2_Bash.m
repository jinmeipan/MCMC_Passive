%Run2_Bash
%check before use it!

model='iba';




%% 'how to name the new r.mat'?  - A: mcmci-empj.mat
%where,i is basic setting plan from 1 to 8
%j is combination of site+month, as:
% site_string={'sdkl','sdkl','sdkl','sdkl','sdkl','churchill','churchill','churchill','churchill','clpx'};
% month_string=[11,12,1,2,3,    1,2,3,4,    2];


planc=1:8;
sitemonc=1:10;

folder='/Users/office/Downloads/MCMC/MCMCRunData_V2/MCMC_Test2freq/';

cd(folder);

% for i=1:length(planc)
for i=2:2
    for j=1:5 %length(sitemonc)

        filename_name=[folder,'MCMC_',num2str(i),'/',model,num2str(j),'/FILENAME.txt'];
    
    r=MCMCRun2;
    r.filename=filename_name;
    
    if(strcmp(model,'hut')==1)
        r.lyrplan=1;
    else
        r.lyrplan=99;%1=use modelSelection to select number of layers; 99=contrained minimum number of layers to 2; 
                     %2=use r.nlyr_choose
    end
    
    
    r.lyrplan=99;
% 	r.nlyr_choose=1; %write when lyplan=2
    
   
    r.run_location=0; %=1: on supper computer,=0 on desktop
    
%     if(strcmp(model,'hut')==1)
%         if( (i==4 | i==5 | i==6 | i==7) & j==3)
%             r.run_location=0;
%         end
%         if(i==2)
%             r.run_location=0;
%         end
%     end
%     
%     if(strcmp(model,'iba')==1)
%         if( (i==2 & j==3) | (i==3 & j==4) |(i==4 & j==8) |(i==4 & j==9))
%             r.run_location=0;
%         end
%     end
%     
%     if(i==3)
%         r.run_location=0;
%     end
    
    r.Max_Lyr=2;
    r=r.Main;
    
    
    name=['mcmc',num2str(i),'-',model,num2str(j)]
    
    
	folder2='/Users/office/Downloads/MCMC/MCMCRunData_V2/';
	save([folder2,'rfile3-test2freq/',name,'.mat'],'r');
    end
end