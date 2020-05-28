clear all

addpath(genpath('impulse_response_materials/Konur'))

tic
%
trlN = 100;
Fs = 1000;  % sample rate
rateN = [5 10 25 50 80 100]/Fs; 
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

%% ======== SPIKE =============
%% extract data for calculation
for kk = 5 %1:length(rateN)
    clearvars -except kk rateN Fs lfp* x* trlN params trln dataLn
    rate = rateN(kk);
    mu=mean(x1,1);
    v=std(x1,0,1);
    x1=(x1-repmat(mu,dataLn,1))./(repmat(v,dataLn,1));
    mu=mean(x2,1);
    v=std(x2,0,1);
    x2=(x2-repmat(mu,dataLn,1))./(repmat(v,dataLn,1));
    mu=mean(x3,1);
    v=std(x3,0,1);
    x3=(x3-repmat(mu,dataLn,1))./(repmat(v,dataLn,1));
    mu=mean(x4,1);
    v=std(x4,0,1);
    x4=(x4-repmat(mu,dataLn,1))./(repmat(v,dataLn,1));
    mu=mean(x5,1);
    v=std(x5,0,1);
    x5=(x5-repmat(mu,dataLn,1))./(repmat(v,dataLn,1));

    %%% generating spikes
    spike1=struct('times',{});
    spike2=struct('times',{});
    spike3=struct('times',{});
    spike4=struct('times',{});
    spike5=struct('times',{});
    for i=1:size(x1,2)
        %%%spike 1
        spike=x1(:,i)>norminv(1-rate);
        s1(i,:)=spike;
        idx_1 =find(spike==1);
        spiketime=idx_1/Fs;
        spike1(i).times=spiketime;
        length1(i)=size(spiketime,1);

        %%%% spike 2
        spike=x2(:,i)>norminv(1-rate);
        s2(i,:)=spike;
        idx_1 =find(spike==1);
        spiketime=idx_1/Fs;
        spike2(i).times=spiketime;
        length2(i)=size(spiketime,1);

        %%%% spike 3
        spike=x3(:,i)>norminv(1-rate);
        s3(i,:)=spike;
        idx_1 =find(spike==1);
        spiketime=idx_1/Fs;
        spike3(i).times=spiketime;
        length3(i)=size(spiketime,1);
        
        %%%% spike 4
        spike=x4(:,i)>norminv(1-rate);
        s4(i,:)=spike;
        idx_1 =find(spike==1);
        spiketime=idx_1/Fs;
        spike4(i).times=spiketime;
        length4(i)=size(spiketime,1);
        
        %%%% spike 3
        spike=x5(:,i)>norminv(1-rate);
        s5(i,:)=spike;
        idx_1 =find(spike==1);
        spiketime=idx_1/Fs;
        spike5(i).times=spiketime;
        length5(i)=size(spiketime,1);
    end
    
    spktr1 = double(s1);
    spktr2 = double(s2);
    spktr3 = double(s3);
    spktr4 = double(s4);
    spktr5 = double(s5);

    %% 1. SPIKE: weighted impulseest - multitaper FFT
    for i = 1:4
        for j = (i+1):5
            eval(['spktri=spktr',num2str(i),'; spktrj=spktr',num2str(j),';']);
            for trl = trln
                [~,~,~,~,P1,P2,f]=coherencypbImp(spktri(trl,:)',spktrj(trl,:)',params); % spike-spike
                datafn1 = iddata(mean(P2,2), mean(P1,2), 1/params.Fs, 'Frequency', f*2*pi);
                datafn2 = iddata(mean(P1,2), mean(P2,2), 1/params.Fs, 'Frequency', f*2*pi);
                try
                sys1fn1 = impulseest(datafn1);       % with TC kernel
                sys1fn2 = impulseest(datafn2);
                
                hf1(:,trl)=sys1fn1.num;
                hf2(:,trl)=sys1fn2.num;
                catch
                end
                clear sys2fn P1 P2
            end
            hf(i,j,:) = mean(hf2,2);
            hf(j,i,:) = mean(hf1,2);

    %% 2. SPIKE: transfer function Txy = Cxy/Cxx

            [C12,phi,P12,P21,P1,P2,ff]=coherencypb(spktri',spktrj',params);
            Sp(i,j,:)= P12; Sp(j,i,:)= P21;
            Sp(i,i,:)= P1;  Sp(j,j,:)= P2;
            % using wilson
            [Snew,Hnew,Znew] = sfactorization_wilsonMG(Sp([i,j],[i,j],:),ff,Fs);
            irpw(i,j,:) = inh_excit1(squeeze(Hnew(1,2,:)));
            irpw(j,i,:) = inh_excit1(squeeze(Hnew(2,1,:)));
            temp = inh_excit(Hnew);
            irpw1(j,i,:) = squeeze(temp(2,1,:));
            irpw1(i,j,:) = squeeze(temp(1,2,:)); clear temp;
            clear Snew Hnew Znew
            % using Txy
            txy1 = P21./P1;
            txy2 = P12./P2;
            bb(i,j,:) = inh_excit1(txy2);
            bb(j,i,:) = inh_excit1(txy1);
        end
    end

    %% 3. SPIKE: spectral matrix factoriaztion

    [Snew1,Hnew1,Znew1] = sfactorization_wilsonMG(Sp,ff,Fs);
    temp = inh_excit(Hnew1);
    
    for i =1:4
        for j=i+1:5
            irpwS(j,i,:) = inh_excit1(squeeze(Hnew1(j,i,:)));
            irpwS(i,j,:) = inh_excit1(squeeze(Hnew1(i,j,:)));
            irpwS1(i,j,:) = squeeze(temp(i,j,:));
            irpwS1(j,i,:) = squeeze(temp(j,i,:)); 
        end
    end
    clear temp;
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
xlim([0 30]);title([num2str(n),'->',num2str(m)]);
line([0 50], [0 0], 'linestyle','--','color','black');
ylabel('amplitude'); %ylim([-0.08 0.02])
hold off

subplot(5,1,2);
hold on
m=3; n=1;
plot(0:size(hf,3)-1,squeeze(hf(m,n,:)),'-ko'); % impulseest
plot(0:size(irpw,3)-1,squeeze(irpw(m,n,:)),'-ro'); % wilson 
plot(0:size(irpw1,3)-1,squeeze(irpw1(m,n,:)),'-mo'); % wilson (VAR) 
plot(0:size(irpwS,3)-1,squeeze(irpwS(m,n,:)),'-go'); % wilson multivar
plot(0:size(irpwS1,3)-1,squeeze(irpwS1(m,n,:)),'-co'); % wilson multivar (VAR)
plot(0:size(bb,3)-1,squeeze(bb(m,n,:)),'-bo'); % txy
xlim([0 30]);title([num2str(n),'->',num2str(m)]);
line([0 50], [0 0], 'linestyle','--','color','black');
ylabel('amplitude');%ylim([-0.08 0.02])
hold off

subplot(5,1,3);
hold on
m=4; n=1;
plot(0:size(hf,3)-1,squeeze(hf(m,n,:)),'-ko'); % impulseest
plot(0:size(irpw,3)-1,squeeze(irpw(m,n,:)),'-ro'); % wilson 
plot(0:size(irpw1,3)-1,squeeze(irpw1(m,n,:)),'-mo'); % wilson (VAR) 
plot(0:size(irpwS,3)-1,squeeze(irpwS(m,n,:)),'-go'); % wilson multivar
plot(0:size(irpwS1,3)-1,squeeze(irpwS1(m,n,:)),'-co'); % wilson multivar (VAR)
plot(0:size(bb,3)-1,squeeze(bb(m,n,:)),'-bo'); % txy
xlim([0 30]);title([num2str(n),'->',num2str(m)]);
line([0 50], [0 0], 'linestyle','--','color','black');
ylabel('amplitude')
%ylim([-0.5 0.5])
hold off

subplot(5,1,4);
hold on
m=5; n=4;
plot(0:size(hf,3)-1,squeeze(hf(m,n,:)),'-ko'); % impulseest
plot(0:size(irpw,3)-1,squeeze(irpw(m,n,:)),'-ro'); % wilson 
plot(0:size(irpw1,3)-1,squeeze(irpw1(m,n,:)),'-mo'); % wilson (VAR) 
plot(0:size(irpwS,3)-1,squeeze(irpwS(m,n,:)),'-go'); % wilson multivar
plot(0:size(irpwS1,3)-1,squeeze(irpwS1(m,n,:)),'-co'); % wilson multivar (VAR)
plot(0:size(bb,3)-1,squeeze(bb(m,n,:)),'-bo'); % txy
xlim([0 30]);title([num2str(n),'->',num2str(m)]);
line([0 50], [0 0], 'linestyle','--','color','black');
ylabel('amplitude')
%ylim([-0.5 0.5])
hold off

subplot(5,1,5);
hold on
m=4; n=5;
plot(0:size(hf,3)-1,squeeze(hf(m,n,:)),'-ko'); % impulseest
plot(0:size(irpw,3)-1,squeeze(irpw(m,n,:)),'-ro'); % wilson 
plot(0:size(irpw1,3)-1,squeeze(irpw1(m,n,:)),'-mo'); % wilson (VAR) 
plot(0:size(irpwS,3)-1,squeeze(irpwS(m,n,:)),'-go'); % wilson multivar
plot(0:size(irpwS1,3)-1,squeeze(irpwS1(m,n,:)),'-co'); % wilson multivar (VAR)
plot(0:size(bb,3)-1,squeeze(bb(m,n,:)),'-bo'); % txy
xlim([0 30]);title([num2str(n),'->',num2str(m)]);
line([0 50], [0 0], 'linestyle','--','color','black');
ylabel('amplitude');xlabel('Time(msec)')
hold off

%%
hold on
legend('Trial-averaged Impulseest of multi-taper FFT)',...
    'Transfer function: Bivariate Wilson Factorization',...
    'Transfer function: Bivariate Wilson Factorization(VAR coeff)',...
    'Transfer function: Multivariate Wilson Factorization',...
    'Transfer function: Multivariate Wilson Factorization(VAR coeff)',...
    'Transfer function: Txy=Pxy/Pxx','Location','northeast');


sgtitle({'5 Node Model';'Spike-Spike';['rate = ',num2str(rate*1000),' Hz']});

pos = get(gcf,'pos');
set(gcf,'pos',[pos(1) pos(2) 8 10])

hold off
toc