function simulation_permutation(data_type,model)

% By Konuralp Bayrak
% Drexel University 2020

fmax =  [250  300  350  400  450  500]; % for lfp
rates=  [5    10   20   50   80   100]; % for spk

if strcmp(data_type,'lfp')
    for maxHz = 1:length(fmax)
        if isfile(['lfp_',model,'_',num2str(fmax(maxHz)),'Hz_thresh.mat'])
            continue
        else
            clearvars -except fmax maxHz model
            eval(['load lfp_',model,'_',num2str(fmax(maxHz)),'Hz.mat']);  
            params.fpass =[0 fmax(maxHz)];
            for nn=1:500
                [C12,phi,P12,P21,P1,P2,f]=coherencyc(lfp1(:,randperm(100)),lfp2,params);
                S1(1,2,:)=P12;
                S1(2,1,:)=P21;
                S1(1,1,:)=P1;S1(2,2,:)=P2;
                [Snew,Hnew,Znew]=sfactorization_wilsonMG(S1,f,Fs);
                h12N(:,nn) = inh_excit1(squeeze(Hnew(2,1,:)));
                h21N(:,nn) = inh_excit1(squeeze(Hnew(1,2,:)));
            end
            eval(['save lfp_',model,'_', num2str(fmax(maxHz)),'Hz_thresh.mat']); 
        end
    end
elseif strcmp(data_type,'spk')
    for mm = 1:length(rates)
        if isfile(['spk_',model,'_',num2str(rates(mm)),'Hz_thresh.mat'])
            continue
        else
            clearvars -except rates mm model
            eval(['load spk_',model,'_',num2str(rates(mm)),'Hz.mat']);  

            for nn=1:500
                [C12,phi,P12,P21,P1,P2,f]=coherencypt(spike1(randperm(100)),spike2,params);
                S1(1,2,:)=P12;
                S1(2,1,:)=P21;
                S1(1,1,:)=P1;S1(2,2,:)=P2;
                [Snew,Hnew,Znew]=sfactorization_wilsonMG(S1,f,Fs);
                h12N(:,nn) = inh_excit1(squeeze(Hnew(2,1,:)));
                h21N(:,nn) = inh_excit1(squeeze(Hnew(1,2,:)));
            end
            eval(['save spk_',model,'_',num2str(rates(mm)),'Hz_thresh.mat']); 
        end
    end
end



