%% READ ME: Set MATLAB path to the current directory to run
clear all
addpath(genpath('base'))

tic

%% model parameters


A=[-0.5 0 0 0 0, 0 0 0 0 0, 0 0 0 0 0;...
    0 0 0 0 0, -0.5 0 0 0 0, 0 0 0 0 0;...
    0 0 0 0 0, 0 0 0 0 0, 0.4 0 0 0 0;...
    0 0 0 -0.5 -0.5, 0.5 0 0 0 0, 0 0 0 0 0;...
    0 0 0 0.5 -0.5, 0 0 0 0 0, 0 0 0 0 0];

%{
A=[-0.5 0 0 0 0, 0 0 0 0 0, 0 0 0 0 0, 0 0 0 0 0, 0 0 0 0 0;...
    0 0 0 0 0, -0.5 0 0 0 0, 0 0 0 0 0, 0 0 0 0 0, 0 0 0 0 0.5;...
    0 0 0 0 0, 0 0 0 0 0, 0 0.5 0 0 0, 0 0 0 0 0, 0 0 0 0 0;...
    0 0 0 0 0, 0 0 0 0 0, 0 0.5 0 0 0, 0 0 0 0 0, 0 0 0 0 0;...
    0 0 0 0 0, 0 0 0 0 0, 0 0 0 0 0, 0 0 -0.5 -0.5 0, 0 0 0 0 0];
%}
%{
A=[-0.95*sqrt(2) 0 0 0 0, 0.9025 0 0 0 0, 0 0 0 0 0;...
    0 0 0 0 0, -0.5 0 0 0 0, 0 0 0 0 0;...
    0 0 0 0 0, 0 0 0 0 0, 0.4 0 0 0 0;...
    0 0 0 -0.25*sqrt(2) -0.25*sqrt(2), 0.5 0 0 0 0, 0 0 0 0 0;...
    0 0 0 0.25*sqrt(2) -0.25*sqrt(2), 0 0 0 0 0, 0 0 0 0 0];
%}

%C=diag([0.6 0.5 0.4 0.4 0.6]);
trlN = 100;

C=eye(5);
[mar] = mar_init(A, C);
dataLn = 1000;

for i=1:100
    temp=mar_gen (mar,dataLn);
    x1(:,i)=temp(1,:)'; 
    x2(:,i)=temp(2,:)';
    x3(:,i)=temp(3,:)';
    x4(:,i)=temp(4,:)';
    x5(:,i)=temp(5,:)';
end
lfp1=x1; lfp2=x2; lfp3=x3; lfp4=x4; lfp5=x5;

%% extract data for calculation

% lfp1=lfp1(101:100+dataLn,1:trlN);
% lfp2=lfp2(101:100+dataLn,1:trlN);
% lfp3=lfp3(101:100+dataLn,1:trlN);
% x1=x1(101:100+dataLn,1:trlN);
% x2=x2(101:100+dataLn,1:trlN);
% x3=x3(101:100+dataLn,1:trlN);
trln = 1:100;     % selected trial

%% params settings
Fs=1000;
params.Fs= Fs; % sampling frequency
params.tapers=[5 9]; % taper parameters
params.trialave=1;
Fs = params.Fs;
params.fpass=[0 Fs/2];
%% 1. LFP: impulse - Wilson factorization
for i =1:4
    for j=(i+1):5
        
        eval(['lfpi=lfp',num2str(i),'; lfpj=lfp',num2str(j),';']);
        for trl = 1:100
            data1 = iddata(lfpi(:,trl), lfpj(:,trl), 1/params.Fs);  % output, input
            data2 = iddata(lfpj(:,trl), lfpi(:,trl), 1/params.Fs);
            opt = impulseestOptions('RegulKernel', 'none'); 
            sys1 = impulseest(data1, opt);  % without kernel
            sys2 = impulseest(data2, opt);  % without kernel
            h1(:,trl) = sys1.num;
            h2(:,trl) = sys2.num;
        end
        hf(i,j,:) = mean(h1,2);
        hf(j,i,:) = mean(h2,2);
        
		[C,phi,P12,P21,P1,P2,ff]= coherencyc(lfpi,lfpj,params);
		Sp(i,j,:)= P12;
        Sp(j,i,:)= P21;
        Sp(i,i,:)= P1;
        Sp(j,j,:)= P2;
        [Snew,Hnew,Znew] = sfactorization_wilsonMG(Sp([i,j],[i,j],:),ff,Fs);
        irpw(i,j,:) = inh_excit1(squeeze(Hnew(1,2,:)));
        irpw(j,i,:) = inh_excit1(squeeze(Hnew(2,1,:)));
        temp = inh_excit(Hnew);
        irpw1(i,j,:) = squeeze(temp(1,2,:));
        irpw1(j,i,:) = squeeze(temp(2,1,:)); clear temp;
        % using Txy
        txy1 = P21./P1;
        txy2 = P12./P2;
        bb(i,j,:) = inh_excit1(txy2);
        bb(j,i,:) = inh_excit1(txy1);
        clear Snew Hnew Znew txy1 txy2 
    end
end

%% pseudo pairwise using multivariate factorization

[Snew1,Hnew1,Znew1] = sfactorization_wilsonMG(Sp,ff,Fs);
temp = inh_excit(Hnew1);

for i =1:4
    for j =(i+1):5
        irpwS(i,j,:) = inh_excit1(squeeze(Hnew1(i,j,:)));
        irpwS(j,i,:) = inh_excit1(squeeze(Hnew1(j,i,:)));
        irpwS1(i,j,:) = squeeze(temp(i,j,:));
        irpwS1(j,i,:) = squeeze(temp(j,i,:)); 
    end
end

        
%% generate plots for report 
figure('units','inches')
subplot(5,1,1);
hold on
m=2; n=1;
plot(0:size(hf,3)-1,squeeze(hf(m,n,:)),'-ko'); % impulseest
plot(0:size(irpw,3)-1,squeeze(irpw(m,n,:)),'-ro'); % wilson 
plot(0:size(irpw1,3)-1,squeeze(irpw1(m,n,:)),'-mo'); % wilson (VAR) 
plot(0:size(irpwS,3)-1,squeeze(irpwS(m,n,:)),'-go'); % wilson multivar
plot(0:size(irpwS1,3)-1,squeeze(irpwS1(m,n,:)),'-co'); % wilson multivar (VAR)
plot(0:size(bb,3)-1,squeeze(bb(m,n,:)),'-bo'); % txy
%line([0 50], [0 0], 'linestyle','--','color','black');
xlim([0 30]);title([num2str(n),'->',num2str(m)]);
ylabel('amplitude')
hold off

subplot(5,1,2);
hold on
m=3; n=1;
plot(0:size(hf,3)-1,squeeze(hf(m,n,:)),'-ko'); % impulse 
plot(0:size(irpw,3)-1,squeeze(irpw(m,n,:)),'-ro'); % wilson 
plot(0:size(irpw1,3)-1,squeeze(irpw1(m,n,:)),'-mo'); % wilson (VAR) 
plot(0:size(irpwS,3)-1,squeeze(irpwS(m,n,:)),'-go'); % wilson multivar
plot(0:size(irpwS1,3)-1,squeeze(irpwS1(m,n,:)),'-co'); % wilson multivar (VAR)
plot(0:size(bb,3)-1,squeeze(bb(m,n,:)),'-bo'); % txy
%line([0 50], [0 0], 'linestyle','--','color','black');
xlim([0 30]);title([num2str(n),'->',num2str(m)]);
ylabel('amplitude')
hold off

subplot(5,1,3);
hold on
m=4; n=1;
plot(0:size(hf,3)-1,squeeze(hf(m,n,:)),'-ko'); % impulse
plot(0:size(irpw,3)-1,squeeze(irpw(m,n,:)),'-ro'); % wilson 
plot(0:size(irpw1,3)-1,squeeze(irpw1(m,n,:)),'-mo'); % wilson (VAR) 
plot(0:size(irpwS,3)-1,squeeze(irpwS(m,n,:)),'-go'); % wilson multivar
plot(0:size(irpwS1,3)-1,squeeze(irpwS1(m,n,:)),'-co'); % wilson multivar (VAR)
plot(0:size(bb,3)-1,squeeze(bb(m,n,:)),'-bo'); % txy
%line([0 50], [0 0], 'linestyle','--','color','black');
xlim([0 30]);title([num2str(n),'->',num2str(m)]);
ylabel('amplitude')
hold off

subplot(5,1,4);
hold on
m=5; n=4;
plot(0:size(hf,3)-1,squeeze(hf(m,n,:)),'-ko'); % impulse
plot(0:size(irpw,3)-1,squeeze(irpw(m,n,:)),'-ro'); % wilson 
plot(0:size(irpw1,3)-1,squeeze(irpw1(m,n,:)),'-mo'); % wilson (VAR) 
plot(0:size(irpwS,3)-1,squeeze(irpwS(m,n,:)),'-go'); % wilson multivar
plot(0:size(irpwS1,3)-1,squeeze(irpwS1(m,n,:)),'-co'); % wilson multivar (VAR)
plot(0:size(bb,3)-1,squeeze(bb(m,n,:)),'-bo'); % txy
%line([0 50], [0 0], 'linestyle','--','color','black');
xlim([0 30]);title([num2str(n),'->',num2str(m)]);
ylabel('amplitude')
hold off

subplot(5,1,5);
hold on
m=4; n=5;
plot(0:size(hf,3)-1,squeeze(hf(m,n,:)),'-ko'); % impulse 
plot(0:size(irpw,3)-1,squeeze(irpw(m,n,:)),'-ro'); % wilson 
plot(0:size(irpw1,3)-1,squeeze(irpw1(m,n,:)),'-mo'); % wilson (VAR) 
plot(0:size(irpwS,3)-1,squeeze(irpwS(m,n,:)),'-go'); % wilson multivar
plot(0:size(irpwS1,3)-1,squeeze(irpwS1(m,n,:)),'-co'); % wilson multivar (VAR)
plot(0:size(bb,3)-1,squeeze(bb(m,n,:)),'-bo'); % txy
%line([0 50], [0 0], 'linestyle','--','color','black');
xlim([0 30]);title([num2str(n),'->',num2str(m)]);
ylabel('amplitude');xlabel('Time(msec)')
hold off

%%
hold on
legend('Trial-averaged Impulseest of LFP mttFFT)',...
    'Transfer function: Bivariate Wilson Factorization',...
    'Transfer function: Bivariate Wilson Factorization(VAR coeff)',...
    'Transfer function: Multivariate Wilson Factorization',...
    'Transfer function: Multivariate Wilson Factorization(VAR coeff)',...
    'Transfer function: Txy=Pxy/Pxx','Location','northeast');


sgtitle({'5 Node Model';'LFP-LFP';['0-',num2str(Fs/2),' Hz']});

pos = get(gcf,'pos');
set(gcf,'pos',[pos(1) pos(2) 8 10])

hold off

toc

