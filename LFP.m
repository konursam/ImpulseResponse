%% READ ME: Set MATLAB path to the current directory to run

addpath(genpath('base'))
%%  only using rate 20Hz
%%% (1)Generating 3 channels of LFP, spikes
%%% generating continuous simulated data based on AR model (2 channels)
clear all
tic
%% model parameters
A=[-0.5 0 0 0 ; 0 -0.5 -0.7 0]; % excitatory
%A=[-0.5 0 0 0 0 0;0 -0.5 0 0 0.7 0]; % inhibitory
%A= -[0.5 0 0 0 0 0; 0.5 0.5 0 0 -0.9 0]; % mixed 1
%A= -[0.5 0 0 0 0 0; -0.5 0.5 0 0 0.9 0]; % mixed 2
%A= -[0.5 0 0, 0 0 0; 0.5 0 0, -0.5 0 0]; % mixed equal kb

%% 
trlN = 100; % number of trials
%m=3;C = diag([noi(m) 1]); % change lfp1 variance
%C = diag([1 noi(m)]); % change lfp2 variance
C = eye(2);
[mar] = mar_init(A, C); 
lyap(mar,1000);
%dataL = 3000;
dataL = 5000;
clear x1 x2 lfp1 lfp2
for i=1:trlN % trlN number of trials
    temp=mar_gen (mar, dataL); % Generates dataL sampes of time series from MAR model
    x1(:,i)=temp(1,:)';
    x2(:,i)=temp(2,:)';
end   
lfp1=x1;lfp2=x2;

%% extract data for calculation
dataLn = 1000; % data length
lfp1=   lfp1(101:100+dataLn,1:trlN);
lfp2=   lfp2(101:100+dataLn,1:trlN);
x1=     x1(101:100+dataLn,1:trlN);
x2=     x2(101:100+dataLn,1:trlN);

%% parameters for multitaper spectrum
params.Fs=1000; % sampling frequency
%params.tapers=[4 7]; % taper parameters
params.tapers=[4 7]; % taper parameters
params.trialave=1;
params.pad = 0;
Fs = params.Fs;
Ts = 1/Fs;
%% Iterate at different FMAX


fmax = [Fs/2-250 Fs/2-200 Fs/2-150 Fs/2-100 Fs/2-50 Fs/2];
figure('units','inches')
for maxHz = 1:length(fmax)
    params.fpass =[0 fmax(maxHz)];

    %% using MATLAB impulseeest
    % ============= LFP================
    %1. time domain (standard)
    for trl = 1:trlN
        data = iddata(lfp2(:,trl), lfp1(:,trl), 1/params.Fs);  % output, input
        
        opt = impulseestOptions('RegulKernel', 'none'); 
        sys1 = impulseest(data, opt);  % without kernel
        h1(:,trl) = sys1.num; 
        
        sys2 = impulseest(data);       % with TC kernel
        h2(:,trl) = sys2.num;
    end

    %% full frequency range 
    clear S
    [C12,phi,P12,P21,P1,P2,f]=coherencyc(lfp1,lfp2,params);
    S(1,2,:)=P12;
    S(2,1,:)=P21;
    S(1,1,:)=P1;
    S(2,2,:)=P2;

    [Snew,Hnew,Znew]=sfactorization_wilsonMG(S,f,Fs);
    temp = inh_excit(Hnew); % Matrix inversion of Hnew
    ir12s1= squeeze(temp(2,1,:));clear temp
    [ir12s] = inh_excit1(squeeze(Hnew(2,1,:)));

    %%%%%% Wilson factorization with full 0-500hz frequency range
    %{
    params.fpass =[0 500];
    clear S
    [C12,phi,P12,P21,P1,P2,f]=coherencyc(lfp1,lfp2,params);
    S(1,2,:)=P12;
    S(2,1,:)=P21;
    S(1,1,:)=P1;
    S(2,2,:)=P2;

    [Snew,Hnew,Znew]=sfactorization_wilsonMG(S,f,Fs);
    [ir12s2] = inh_excit1(squeeze(Hnew(2,1,:)));

    %Txy = squeeze(Snew(1,2,:))./squeeze(Snew(1,1,:));
    %Txy12 = inh_excit1(Txy); 
    %}
    %title(['0-',num2str(fff),' Hz']);

    %%% LFP: Txy
    txy = P21./P1;
    %[blfp,a] = invfreqz(txy,f*2*pi/Fs,70,0);
    blfp = inh_excit1(txy);
    
    if A == [-0.5 0 0 0 ; 0 -0.5 -0.7 0]
        eval(['save lfp_ext_', num2str(fmax(maxHz)),'Hz.mat;'])
    elseif A == [-0.5 0 0 0 0 0;0 -0.5 0 0 0.7 0]
        eval(['save lfp_inh_', num2str(fmax(maxHz)),'Hz.mat;'])
    elseif A == -[0.5 0 0, 0 0 0; 0.5 0 0, -0.5 0 0]
        eval(['save lfp_mixed_equal_', num2str(fmax(maxHz)),'Hz.mat;'])
    elseif A == -[0.5 0 0 0 0 0; 0.5 0.5 0 0 -0.9 0]
        eval(['save lfp_mixed1_', num2str(fmax(maxHz)),'Hz.mat;'])
    elseif A == -[0.5 0 0 0 0 0; -0.5 0.5 0 0 0.9 0]
        eval(['save lfp_mixed2_', num2str(fmax(maxHz)),'Hz.mat;'])
    end
    %end

    %% plot


    subplot(round(length(fmax)/2),2,maxHz)
    hold on
    plot(0:size(h1,1)-1, mean(h1,2), '-k.'); % impulseest
    plot(0:length(ir12s)-1,ir12s,'-r.'); %  wilson
    plot(0:length(ir12s1)-1,ir12s1,'-m.'); % wilson VAR coef
    %plot(0:length(ir12s2)-1,ir12s2,'-g.');
    plot(0:length(blfp)-1,blfp,'-b.'); % Txy
    line([0 50], [0 0], 'linestyle','--','color','black');
    xlabel('Time(msec)'); 
    ylabel('Amplitude');
    
    xlim([0 30]); %ylim([-0.2 1]);
    title({['0-',num2str(fmax(maxHz)),' Hz']},'FontWeight','Normal');
    
end


    
if A == [-0.5 0 0 0 ; 0 -0.5 -0.7 0]
    sgtitle({'Excitation Model';'LFP-LFP'});
elseif A == [-0.5 0 0 0 0 0;0 -0.5 0 0 0.7 0]
    sgtitle({'Inhibition Model';'LFP-LFP'});
elseif A == -[0.5 0 0, 0 0 0; 0.5 0 0, -0.5 0 0]
    sgtitle({'Mixed Equal Model';'LFP-LFP'});
elseif A == -[0.5 0 0 0 0 0; 0.5 0.5 0 0 -0.9 0]
    sgtitle({'Mixed Model 1';'LFP-LFP'});
elseif A == -[0.5 0 0 0 0 0; -0.5 0.5 0 0 0.9 0]
    sgtitle({'Mixed Model 2';'LFP-LFP'});
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


%% Focused plot on Fmax = Fs/2

figure
%subplot(211)
hold on

plot(0:size(h1,1)-1, mean(h1,2), '-ko'); % impulseest trial ave
plot(0:length(ir12s)-1,ir12s,'-ro'); % wilson
%area(0:length(ir12s)-1,ir12s);
plot(0:length(ir12s1)-1,ir12s1,'-mo'); % wilson VAR coeff
plot(0:length(blfp)-1,blfp,'-bo'); % Txy

line([0 50], [0 0], 'linestyle','--','color','black');
xlabel('Time(sec)'); 
ylabel('Amplitude');
legend({'Trial-averaged impulseest of multi-taper FFT',...
'Transfer function: Wilson Factorization',...
'Transfer function: Wilson Factorization (VAR coeff)',...
'Transfer function: Txy=Pxy/Pxx'},'Location','best');
xlim([0 30]);% ylim([-1 1]);
Q = trapz(0:length(ir12s)-1,ir12s);
%title({['Area = ',num2str(Q)]},'FontWeight','Normal');

if A == [-0.5 0 0 0 ; 0 -0.5 -0.7 0]
    title({'Excitation Model';'LFP-LFP';['0-',num2str(fmax(maxHz)),' Hz']});
elseif A == [-0.5 0 0 0 0 0;0 -0.5 0 0 0.7 0]
    title({'Inhibition Model';'LFP-LFP';['0-',num2str(fmax(maxHz)),' Hz']});
elseif A == -[0.5 0 0, 0 0 0; 0.5 0 0, -0.5 0 0]
    title({'Mixed Equal Model';'LFP-LFP';['0-',num2str(fmax(maxHz)),' Hz']});
elseif A == -[0.5 0 0 0 0 0; 0.5 0.5 0 0 -0.9 0]
    title({'Mixed Model 1';'LFP-LFP';['0-',num2str(fmax(maxHz)),' Hz']});
elseif A == -[0.5 0 0 0 0 0; -0.5 0.5 0 0 0.9 0]
    title({'Mixed Model 2';'LFP-LFP';['0-',num2str(fmax(maxHz)),' Hz']});
end

hold off
toc