
clear all
fmax = [250 300 350 400 450 500];

%%
% __________________
% Choose your model:
% ------------------

%model = 'ext';          mod_full = 'Excitation Model';
model = 'inh';         mod_full = 'Inhibition Model';
%model = 'mixed1';      mod_full = 'Mixed Model 1';
%model = 'mixed2';      mod_full = 'Mixed Model 2';
%model = 'mixed_equal'; mod_full = 'Mixed Equal Model';

%% Permutation generator

simulation_permutation('lfp',model);

%% Plot Generator

figure;
for maxHz = 1:length(fmax)
    eval(['load lfp_',model,'_',num2str(fmax(maxHz)),'Hz_thresh.mat']);
    
    %h12Nnew = [ir12s,h12N];
    %p = kruskalwallis(h12Nnew,[],'off');
    for i = 1:length(ir12s)
        x = h12N(i,:)';
        lower(i) = norminv(0.05,mean(x),std(x));
        upper(i) = norminv(1-0.05,mean(x),std(x));
    end
    [M,I] = max(abs(ir12s));
    if abs(min(ir12s)) == M
        M = -1*M;
    end
    x = h12N(I,:)';
    pd = fitdist(x, 'Normal');
    
    if sign(M) == 1
        index = linspace(min(x)+min(x)*0.5, max(x)+max(x)*0.1, 1000);
    else
        index = linspace(min(x)+min(x)*0.5, max(x)+max(x)*0.1, 1000);
    end
    
    h(maxHz) = subplot(round(length(fmax)/2),2,maxHz);
    
    hold on
    y = pdf(pd, index);
    plot(index,y,'-r','LineWidth',2);
    line([M M],[0 max(y)+max(y)/2],'color','green','LineWidth',2);
    line([lower(I) lower(I)],[0 max(y)-max(y)/3],'color','blue','LineWidth',1);
    line([upper(I) upper(I)],[0 max(y)-max(y)/3],'color','blue','LineWidth',1);
    title(['fpass = [0 ',num2str(fmax(maxHz)),'] Hz'],...%,['p=',num2str(p)]},...
        'FontWeight','Normal');
    p(maxHz,:) = get(h(maxHz),'position');
    hold off
end

hold on

legend({'Normally Dist. Permutations',...
    'Wilson decomp. IR',...
    '95% confidence int.'},...
    'Location','best');
legend('boxoff')

height = p(1,2) + p(1,4) - p(maxHz,2);
width = p(maxHz,1) + p(maxHz,3) - p(1,1);
h(maxHz+1) = axes('position',[p(maxHz-1,1) p(maxHz-1,2) width height],'visible','off');
h(maxHz+1).XLabel.Visible = 'on';
h(maxHz+1).YLabel.Visible = 'on';

axes(h(maxHz+1))
xlabel('Amplitude');
ylabel({'Pobability Density';'at IR peak'});


sgtitle({[mod_full,' vs Permutations (n=500)'];'LFP-LFP'})


hold off

clearvars -except fmax model mod_full

figure
for maxHz = 1:length(fmax)  
    eval(['load lfp_',model,'_',num2str(fmax(maxHz)),'Hz_thresh.mat']);
    
    for i = 1:length(ir12s)
        x = h12N(i,:)';
        lower(i) = norminv(0.05,mean(x),std(x));
        upper(i) = norminv(1-0.05,mean(x),std(x));
    end
    
    h(maxHz) = subplot(round(length(fmax)/2),2,maxHz);
    hold on
    p1 = plot(0:length(h12N)-1,h12N,'-r.');
    p2 = plot(0:length(ir12s)-1,ir12s,'-g.','LineWidth',1);
    p3 = plot(0:length(lower)-1,lower,'-b','LineWidth',1);
    p4 = plot(0:length(upper)-1,upper,'-b','LineWidth',1);
    xlim([0 30]);
    line([0 30], [0 0], 'linestyle','--','color','black');
    title(['fpass = [0 ',num2str(fmax(maxHz)),'] Hz'],...%['p=',num2str(p)]},...
        'FontWeight','Normal');
    p(maxHz,:) = get(h(maxHz),'position');
    hold off
    %{
    h12Nnew1=mean(h12N,2);
    t = ttest2(ir12s(1:30),h12Nnew1(1:30));
    if t==0
        fprintf('Mixed Model 2 paired t-test for fpass [0 %u ] Hz is not significant\n',fmax(maxHz));
    else
        fprintf('Mixed Model 2 paired t-test for fpass [0 %u ] Hz is significant\n',fmax(maxHz));
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

height = p(1,2) + p(1,4) - p(maxHz,2);
width = p(maxHz,1) + p(maxHz,3) - p(1,1);
h(maxHz+1) = axes('position',[p(maxHz-1,1) p(maxHz-1,2) width height],'visible','off');
h(maxHz+1).XLabel.Visible = 'on';
h(maxHz+1).YLabel.Visible = 'on';

axes(h(maxHz+1))
xlabel('time(msec)');
ylabel('Amplitude');


sgtitle({[mod_full,' vs Permutations (n=500)'];'LFP-LFP'})


hold off
