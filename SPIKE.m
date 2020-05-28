%% READ ME: Set MATLAB path to the current directory to run

addpath(genpath('base'))

clear all
tic
trlN = 100; % number of trials
noi = [0.2 0.5 1 1.5 2 5 100];
cols = ['g','m','k','r','b','m','c'];


%% model parameters
A=[-0.5 0 0 0 ; 0 -0.5 -0.7 0]; % excitatory
%A=[-0.5 0 0 0 0 0;0 -0.5 0 0 0.7 0]; % inhibitory
%A= -[0.5 0 0 0 0 0; 0.5 0.5 0 0 -0.9 0]; % mixed model 1
%A= -[0.5 0 0 0 0 0; -0.5 0.5 0 0 0.9 0]; % mixed model 2 
%A= -[0.5 0 0, 0 0 0; 0.5 0 0, -0.5 -0.5 0]; % mixed equal kb
%for m = 1:7
    %clear x1 x2 lfp1 lfp2
    
%% 
%m=3;C = diag([noi(m) 1]); % change lfp1 variance
%C = diag([1 noi(m)]); % change lfp2 variance
C = diag([1 1]);
[mar] = mar_init(A, C); 
lyap(mar,1000);
%dataL = 3000;
dataL = 5000;
clear x1 x2 lfp1 lfp2
for i=1:trlN 
    temp=mar_gen (mar, dataL);
    x1(:,i)=temp(1,:)';
    x2(:,i)=temp(2,:)';
end   
lfp1=x1;lfp2=x2;

%% extract data for calculation
dataLn = 1000; % data length
lfp1=lfp1(101:100+dataLn,1:trlN);
lfp2=lfp2(101:100+dataLn,1:trlN);
x1=x1(101:100+dataLn,1:trlN);
x2=x2(101:100+dataLn,1:trlN);

%% parameters for multitaper spectrum
params.Fs=1000; % sampling frequency
params.tapers=[4 7]; % taper parameters
params.trialave=1;
params.pad = 0;
Fs = params.Fs;

%rateN= [5 10 20 50 80 100]/Fs;   % different rates 
rateN= [5 10 20 50 80 100]/1000;   % different rates 

params.fpass=[0 Fs/2];


%% =============== spike ===========================
%% iteration through different spike rate

figure('units','inches')
for kk = 1:length(rateN) %5

    clearvars -except kk rateN Fs lfp1 lfp2 x1 x2 trlN params trln dataLn A spike1 spike2 ...
        s1 s2 spktr1 spktr2
    rate = rateN(kk);
    mu=mean(x1,1);
    v=std(x1,0,1);
    x1=(x1-repmat(mu,dataLn,1))./(repmat(v,dataLn,1)); % Channel 1 lfp normalization
    % (x1(i)-xu(i))/std(x1)
    mu=mean(x2,1);
    v=std(x2,0,1);
    x2=(x2-repmat(mu,dataLn,1))./(repmat(v,dataLn,1)); % Channel 2 lfp normalization

    %%% generating spikes
    spike1=struct('times',{});
    spike2=struct('times',{});

    for i=1:size(x1,2) % size(x1,2)=100
        %%% generating spike 1 from x1 & firing rate
        spike=x1(:,i)>norminv(1-rate);
        s1(:,i)=spike;
        idx_1 =find(spike==1);
        spiketime=idx_1/Fs;
        spike1(i).times=spiketime;
        length1(i)=size(spiketime,1);
        %%%% generating spike 2 from x2 & firing rate
        spike=x2(:,i)>norminv(1-rate);
        s2(:,i)=spike;
        idx_1 =find(spike==1);
        spiketime=idx_1/Fs;
        spike2(i).times=spiketime;
        length2(i)=size(spiketime,1);
    end

    spktr1 = double(s1); 
    spktr2 = double(s2);

    %% 1. SPIKE: weighted impulseest - multitaper FFT
    %%% single-trial
    
    for trl = 1:trlN
        [~,~,~,~,P1,P2,f]=coherencypbImp(spktr1(:,trl),spktr2(:,trl),params); % spike-spike
        datafn = iddata(mean(P2,2), mean(P1,2), 1/params.Fs, 'Frequency', f*2*pi);
        try
        sys1fn = impulseest(datafn);       % with TC kernel
        hf2(:,trl)=sys1fn.num;
        catch
        end
    clear sys2fn P1 P2
    end
%}


    %%% (2)Trial average 
    %{
    sumHCov1 = zeros(70,1); sumCov1 = zeros(70,70);
    sumHCov2 = zeros(70,1); sumCov2 = zeros(70,70);
    [C12,phi,P12,P21,P1,P2,f]=coherencypbImp(spktr1,spktr2,params);

    for m = 1 : size(P1,2)
        dataf = iddata(squeeze(mean(P2(:,m,:),3)), squeeze(mean(P1(:,m,:),3)), 1/params.Fs, 'Frequency', f*2*pi);

        opt = impulseestOptions('RegulKernel', 'none');
        sys1f = impulseest(dataf, opt);       % i.without kernel;
        hf1(:,m)=sys1f.num;
        tempCov1 = getcov(sys1f);
        sumHCov1 = sumHCov1 + inv(tempCov1(1:70,1:70))*sys1f.num';
        sumCov1 = sumCov1 + inv(tempCov1(1:70,1:70));


        sys2f = impulseest(dataf);       % ii.with TC kernel;
        hf2(:,m) = sys2f.num;           % smooth
        tempCov2 = getcov(sys2f);
        sumHCov2 = sumHCov2 + inv(tempCov2(1:70,1:70))*sys2f.num';
        sumCov2 = sumCov2 + inv(tempCov2(1:70,1:70));
    end
    wh1 = sumCov1\sumHCov1;  % weighted without kernel
    wh2 = sumCov2\sumHCov2;  % weighted with kernel
%}


    %%% (3)Tapers average
    %{
    sumHCov1 = zeros(70,1); sumCov1 = zeros(70,70);
    sumHCov2 = zeros(70,1); sumCov2 = zeros(70,70);
    [C12,phi,P12,P21,P1,P2,f]=coherencypbImp(spktr1,spktr2,params);

    for m = 1 : size(P1,3)
        dataf = iddata(squeeze(mean(P2(:,:,m),2)), squeeze(mean(P1(:,:,m),2)), 1/params.Fs, 'Frequency', f*2*pi);

        opt = impulseestOptions('RegulKernel', 'none');
        sys1f = impulseest(dataf, opt);       % i.without kernel;
        hf1(:,m)=sys1f.num;
        tempCov1 = getcov(sys1f);
        sumHCov1 = sumHCov1 + tempCov1(1:70,1:70)\(sys1f.num');
        sumCov1 = sumCov1 + inv(tempCov1(1:70,1:70));

        sys2f = impulseest(dataf);       % ii.with TC kernel;
        hf2(:,m) = sys2f.num;           % smooth
        tempCov2 = getcov(sys2f);
        sumHCov2 = sumHCov2 + tempCov2(1:70,1:70)\(sys2f.num');
        sumCov2 = sumCov2 + inv(tempCov2(1:70,1:70));
    end
    wh1 = sumCov1\sumHCov1;
    wh2 = sumCov2\sumHCov2;

%}

    %% 2. SPIKE: transfer function Txy = Cxy/Cxx
    clear S1 txy2
    %[C12,phi,P12,P21,P1,P2,f]=coherencypb(spktr1,spktr2,params);
    [C12,phi,P12,P21,P1,P2,f]=coherencypt(spike1,spike2,params);
    S1(1,2,:)=P12;
    S1(2,1,:)=P21;
    S1(1,1,:)=P1;
    S1(2,2,:)=P2;
    txy2 = P12./P1;
    %[bb,a] = invfreqz(txy2,f*2*pi/Fs,70,0);
    [bb] = inh_excit1(txy2);

    %% 3. SPIKE: spectral matrix factoriaztion
    clear Snew Hnew Znew 
    [Snew,Hnew,Znew]=sfactorization_wilsonMG(S1,f,Fs);
    [Fxy,Fyx,pp1,pp2]= pwcausalr_c(Snew,Hnew,Znew,f,params.Fs);
    h12 = inh_excit1(squeeze(Hnew(1,2,:)));
    temp = inh_excit(Hnew);
    h121=squeeze(temp(1,2,:)); clear temp;
    
    if A == [-0.5 0 0 0 ; 0 -0.5 -0.7 0]
        eval(['save spk_ext_', num2str(rate*Fs),'Hz.mat;'])
    elseif A == [-0.5 0 0 0 0 0;0 -0.5 0 0 0.7 0]
        eval(['save spk_inh_', num2str(rate*Fs),'Hz.mat;'])
    elseif A == -[0.5 0 0, 0 0 0; 0.5 0 0, -0.5 0 0]
        eval(['save spk_mixed_equal_', num2str(rate*Fs),'Hz.mat;'])
    elseif A == -[0.5 0 0 0 0 0; 0.5 0.5 0 0 -0.9 0]
        eval(['save lfp_mixed1_', num2str(rate*Fs),'Hz.mat;'])
    elseif A == -[0.5 0 0 0 0 0; -0.5 0.5 0 0 0.9 0]
        eval(['save lfp_mixed2_', num2str(rate*Fs),'Hz.mat;'])
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot
    
    subplot(round(length(rateN)/2),2,kk);
    %subplot(211)
    hold on
    plot(0:size(hf2,1)-1, mean(hf2,2), '-k.'); % impulseest
    plot(0:length(h12)-1,h12,'-r.'); % Wilson
    plot(0:length(h121)-1,h121,'-m.'); % Wilson VAR coeff
    plot(0:length(bb)-1,bb,'-b.');  % Txy

    xlabel('Time(msec)'); 
    ylabel('Amplitude');
    xlim([0 30]);%ylim([-0.15 0.4]); 
    line([0 30], [0 0], 'linestyle','--','color','black');
    title({['rate = ',num2str(rate*1000),' Hz']},'FontWeight','Normal');
    %title({'Spike-Spike';['spike rate: ',num2str(rate*1000),' Hz']},'FontWeight','Normal');
end

%%% for large subplots
%{
[ax,h1]=suplabel('Time(sec)');
[ax,h2]=suplabel('Amplitude','y');
[ax,h3]=suplabel('Spike-Spike' ,'t');
%}
if A == [-0.5 0 0 0 ; 0 -0.5 -0.7 0]
    sgtitle({'Excitation Model';'Spike-Spike'});
elseif A == [-0.5 0 0 0 0 0;0 -0.5 0 0 0.7 0]
    sgtitle({'Inhibition Model';'Spike-Spike'});
elseif A == -[0.5 0 0, 0 0 0; 0.5 0 0, -0.5 0 0]
    sgtitle({'Mixed Equal Model';'Spike-Spike'});
elseif A == -[0.5 0 0 0 0 0; 0.5 0.5 0 0 -0.9 0]
    sgtitle({'Mixed Model 1';'Spike-Spike'});
elseif A == -[0.5 0 0 0 0 0; -0.5 0.5 0 0 0.9 0]
    sgtitle({'Mixed Model 2';'Spike-Spike'});
end
legend({'impulseest',...
        'H: Wilson Decomp.',...
        'H: Wilson Decomp. (VAR coeff)',...
        'Txy = Pxy/Pxx'},...
        'Location','best');
legend('boxoff')

pos = get(gcf,'pos');
set(gcf,'pos',[pos(1) pos(2) 10 10])

hold off

%% focused plot of highest spike rate: 100Hz
figure
hold on
plot(0:size(hf2,1)-1, mean(hf2,2), '-ko'); % impulseest
plot(0:length(h12)-1,h12,'-ro'); % Wilson
plot(0:length(h121)-1,h121,'-mo'); % Wilson VAR coeff
plot(0:length(bb)-1,bb,'-bo');  % Txy

xlabel('Time(msec)'); 
ylabel('Amplitude');
xlim([0 30]);%ylim([-0.15 0.4]); 
line([0 30], [0 0], 'linestyle','--','color','black');

if A == [-0.5 0 0 0 ; 0 -0.5 -0.7 0]
    title({'Excitation Model';'Spike-Spike';['rate = ',num2str(rate*1000),' Hz']});
elseif A == [-0.5 0 0 0 0 0;0 -0.5 0 0 0.7 0]
    title({'Inhibition Model';'Spike-Spike';['rate = ',num2str(rate*1000),' Hz']});
elseif A == -[0.5 0 0, 0 0 0; 0.5 0 0, -0.5 0 0]
    title({'Mixed Equal Model';'Spike-Spike';['rate = ',num2str(rate*1000),' Hz']});
elseif A == -[0.5 0 0 0 0 0; 0.5 0.5 0 0 -0.9 0]
    title({'Mixed Model 1';'Spike-Spike';['rate = ',num2str(rate*1000),' Hz']});
elseif A == -[0.5 0 0 0 0 0; -0.5 0.5 0 0 0.9 0]
    title({'Mixed Model 2';'Spike-Spike';['rate = ',num2str(rate*1000),' Hz']});
end
legend({'Trial-averaged impulseest of multi-taper FFT',...
    'Transfer function: Wilson Factorization',...
    'Transfer function: Wilson Factorization (VAR coeff)',...
    'Transfer function: Txy=Pxy/Pxx'},'Location','best');

hold off

toc