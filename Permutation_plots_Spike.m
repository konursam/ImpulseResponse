
clear all
rates= [5 10 20 50 80 100];

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

%% Plot Generator

figure;
for mm = 1:length(rates)
    eval(['load spk_',model,'_',num2str(rates(mm)),'Hz_thresh.mat']);
    
    %h12Nnew = [h12,h12N];
    %p = kruskalwallis(h12Nnew,[],'off');
    for i = 1:length(h12)
        x = h12N(i,:)';
        lower(i) = norminv(0.05,mean(x),std(x));
        upper(i) = norminv(1-0.05,mean(x),std(x));
    end
    [M,I] = max(abs(h12));
    if abs(min(h12)) == M
        M = -1*M;
    end
    x = h12N(I,:)';
    pd = fitdist(x, 'Normal');
    
    if sign(M) == 1
        index = linspace(min(x)+min(x)*0.5, max(x)+max(x)*0.1, 1000);
    else
        index = linspace(min(x)+min(x)*0.5, max(x)+max(x)*0.1, 1000);
    end
    
    h(mm) = subplot(round(length(rates)/2),2,mm);
    
    hold on
    y = pdf(pd, index);
    plot(index, y,'-r','LineWidth',2);
    line([M M],[0 max(y)+max(y)/2],'color','green','LineWidth',2);
    line([lower(I) lower(I)],[0 max(y)-max(y)/3],'color','blue','LineWidth',1);
    line([upper(I) upper(I)],[0 max(y)-max(y)/3],'color','blue','LineWidth',1);
    title(['rate = ',num2str(rates(mm)),' Hz'],...%['p=',num2str(p)]},...
        'FontWeight','Normal');
    p(mm,:) = get(h(mm),'position');
    hold off
end

hold on

legend({'Normally Dist. Permutations',...
    'Wilson decomp. IR',...
    '95% confidence int.'},...
    'Location','best');
legend('boxoff')

height = p(1,2) + p(1,4) - p(mm,2);
width = p(mm,1) + p(mm,3) - p(1,1);
h(mm+1) = axes('position',[p(mm-1,1) p(mm-1,2) width height],'visible','off');
h(mm+1).XLabel.Visible = 'on';
h(mm+1).YLabel.Visible = 'on';

axes(h(mm+1))
xlabel('Amplitude');
ylabel({'Pobability Density';'at IR peak'});


sgtitle({[mod_full,' vs Permutations (n=500)'];'Spike-Spike'})


hold off

clearvars -except rates model mod_full

figure
for mm = 1:length(rates)
    eval(['load spk_',model,'_',num2str(rates(mm)),'Hz_thresh.mat']);
    
    for i = 1:length(h12)
        x = h12N(i,:)';
        lower(i) = norminv(0.05,mean(x),std(x));
        upper(i) = norminv(1-0.05,mean(x),std(x));
    end
    
    
    h(mm) = subplot(round(length(rates)/2),2,mm);
    hold on
    p1 = plot(0:length(h12N)-1,h12N,'-r.');
    p2 = plot(0:length(h12)-1,h12,'-g.','LineWidth',1);
    p3 = plot(0:length(lower)-1,lower,'-b','LineWidth',1);
    p4 = plot(0:length(upper)-1,upper,'-b','LineWidth',1);
    xlim([0 30]);
    line([0 30], [0 0], 'linestyle','--','color','black');
    title(['rate = ',num2str(rates(mm)),' Hz'],...%['p=',num2str(p)]},...
        'FontWeight','Normal');
    p(mm,:) = get(h(mm),'position');
    hold off
    %{
    h12Nnew1=mean(h12N,2);
    t = ttest2(h12(1:30),h12Nnew1(1:30));
    if t==0
        fprintf('Mixed Model 2 paired t-test for %u Hz is not significant\n',rates(mm));
    else
        fprintf('Mixed Model 2 paired t-test for %u Hz is significant\n',rates(mm));
    end
    %}
end

hold on

legend([p1(1) p2 p3],...
    {'Randomized Permutations',...
    'Wilson decomp. IR',...
    '95% confidence int.'},...
    'Location','best');
legend('boxoff')

height = p(1,2) + p(1,4) - p(mm,2);
width = p(mm,1) + p(mm,3) - p(1,1);
h(mm+1) = axes('position',[p(mm-1,1) p(mm-1,2) width height],'visible','off');
h(mm+1).XLabel.Visible = 'on';
h(mm+1).YLabel.Visible = 'on';

axes(h(mm+1))
xlabel('time(msec)');
ylabel('Amplitude');


sgtitle({[mod_full,'vs Permutations (n=500)'];'Spike-Spike'})


hold off