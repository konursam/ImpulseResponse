
clearvars

%%
% __________________
% Choose your model:
% ------------------

model = 'ext';          mod_full = 'Excitation Model';
%model = 'inh';         mod_full = 'Inhibition Model';
%model = 'mixed1';      mod_full = 'Mixed Model 1';
%model = 'mixed2';      mod_full = 'Mixed Model 2';
%model = 'mixed_equal'; mod_full = 'Mixed Equal Model';
%% Permutation generator

simulation_permutation('spk',model);

%%
eval(['load spk_',model,'_80Hz_thresh.mat']);

for i = 1:length(h12)
    x = h12N(i,:)';
    lower(i) = norminv(0.05,mean(x),std(x));
    upper(i) = norminv(1-0.05,mean(x),std(x));
end

%% IR figure
% defines x and y coordinates for analysis
x = 0:30;
IR = h12(1:31)';

dIR = gradient(IR);         % gradient
dIR2 = diff(IR)./diff(x);   % diff

% finds max of the gradient 
[M,I] = max(abs(dIR));
if abs(min(dIR)) == M
    M = -1*M;
end

% finds the max of the diff
[M2,I2] = max(abs(dIR2));
if abs(min(dIR2)) == M2
    M2 = -1*M2;
end

% determines polarity of IR from sign of the max GRADIENT
if sign(M) == 1
    polarity = 'positive';
else
    polarity = 'negative';
end
% determines polarity of IR from sign of the max DIFF
if sign(M2) == 1
    polarity2 = 'positive';
else
    polarity2 = 'negative';
end


figure
h1 = subplot(2,1,1);
hold on
plot(x,IR,'-g.');
plot(x,dIR,'-k');
plot(x(1:end-1),dIR2,'-b');
plot(x(I),M,'r*');
plot(x(I2),M2,'r*');

plot(x,h12N(1:31,:)','-r.');
plot(x,lower(1:31)','-b.');
plot(x,upper(1:31)','-b.');
line([0 30], [0 0], 'linestyle','--','color','black');
title(['Spike rate = ',num2str(rates(mm)),' Hz'],...
    'FontWeight','Normal');
sgtitle({[mod_full,' IR'];'Spike-Spike'});
legend('IR','1st derivative of raw IR(gradient)',...
    '1st derivative of raw IR(diff)',...
    ['Maximum value polarity(gradient): ',polarity],...
    ['Maximum value polarity(diff): ',polarity2],...
    'Location','northeast');
p1 = get(h1,'position');
hold off


h2 = subplot(2,1,2);
hold on
plot(x,IR,'-g.');
plot(x,lower(1:31)','-b.');
plot(x,upper(1:31)','-b.');
line([0 30], [0 0], 'linestyle','--','color','black');
title(['Spike rate = ',num2str(rates(mm)),' Hz'],...
    'FontWeight','Normal');
legend('IR','upper 95% CI','lower 95% CI',...
    'Location','northeast');
xlabel('time(msec)');
p2 = get(h2,'position');

height = p1(2) + p1(4) - p2(2);
width = p2(1) + p2(3) - p1(1);
h3 = axes('position',[p2(1) p2(2) width height],'visible','off');
h3.YLabel.Visible = 'on';
axes(h3)
ylabel('Amplitude');
hold off

%% Point test V1
%
% checks each point and its predecesor point to see if within the
% confidence interval. if both points are within the interval, the for-loop
% breaks.

%{
new_IR = zeros(1,length(IR));
new_IR(1:2) = IR(1:2);
for i = 3:length(IR)
    if sign(IR(i)) == 1
        if IR(i) < upper(i) && IR(i-1) < upper(i-1)
            break
        elseif IR(i) > upper(i) || IR(i-1) > upper(i-1)
            new_IR(i) = IR(i);
        end
    else
        if IR(i) > lower(i) && IR(i-1) > lower(i-1)
            break
        elseif IR(i) < lower(i) || IR(i-1) < lower(i-1)
            new_IR(i) = IR(i);
        end
    end
end

areaIR = trapz(x,IR);
area_new_IR = trapz(x,new_IR);

if sign(areaIR) == 1
    polarity_areaIR = 'positive';
else
    polarity_areaIR = 'negative';
end

if sign(area_new_IR) == 1
    polarity_areaIRnew = 'positive';
else
    polarity_areaIRnew = 'negative';
end

figure
subplot(2,1,1);
hold on
plot(x,IR,'-r*');
plot(x,lower(1:31)','-b.');
plot(x,upper(1:31)','-b.');
area(x,IR);
line([0 30], [0 0], 'linestyle','--','color','black');
title(['Spike rate = ',num2str(rates(mm)),' Hz'],...
    'FontWeight','Normal');
legend('IR','upper 95% CI','lower 95% CI',['area = ',num2str(areaIR)],...
    ['polarity est from area : ',polarity_areaIR],...
    'Location','northeast');
hold off

subplot(2,1,2);
hold on
plot(x,new_IR,'-r*');
plot(x,lower(1:31)','-b.');
plot(x,upper(1:31)','-b.');
area(x,new_IR);
line([0 30], [0 0], 'linestyle','--','color','black');
legend('IR','upper 95% CI','lower 95% CI',...
    ['area = ',num2str(area_new_IR)],...
    ['polarity est from area : ',polarity_areaIRnew],...
    'Location','northeast');
hold off

%}

