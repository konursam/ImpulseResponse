
clear all

fmax = [250 300 350 400 450 500];

%% Points to observe

maxpoints = 20;  % choose the number of points you'd like to view
maxpoints = 2*round(maxpoints/2); % makes maxpoints even number

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

simulation_permutation('lfp',model);
%% Plot generator

for maxHz = 1:length(fmax)
    clearvars -except fmax maxHz k maxpoints model mod_full
    eval(['load lfp_',model,'_',num2str(fmax(maxHz)),'Hz_thresh.mat']);
    %h12Nnew = [ir12s,h12N];
    %p = kruskalwallis(h12Nnew,[],'off');
    for i = 1:length(ir12s)
        x = h12N(i,:)';
        lower(i) = norminv(0.05,mean(x),std(x));
        upper(i) = norminv(1-0.05,mean(x),std(x));
    end
    
    figure('units','inches')
    
    subplot((maxpoints/2)+1,2,[1,2]);
    hold on
    p1 = plot(0:length(h12N)-1,h12N,'-r.');
    p2 = plot(0:length(ir12s)-1,ir12s,'-go','LineWidth',2);
    p3 = plot(0:length(lower)-1,lower,'-b','LineWidth',1);
    p4 = plot(0:length(upper)-1,upper,'-b','LineWidth',1);
    xlim([0 30]);
    xlabel('time(msec)');
    ylabel('Amplitude');
    line([0 30], [0 0], 'linestyle','--','color','black');
    
    legend([p1(1) p2 p3],...
        {'Randomized Permutations',...
        'Wilson decomp. IR',...
        '95% confidence int.'},...
        'Location','best');
    legend('boxoff')

    
    for k =2:maxpoints+1
        
        M = ir12s(k);
        
        x = h12N(k,:)';
        pd = fitdist(x, 'Normal');
        index = linspace(min(x), max(x), 1000);

        %figure(1)
        h(k+1) = subplot((maxpoints/2)+1,2,k+1);
        y = pdf(pd, index);
        plot(index, y,'-r','LineWidth',2);
        line([M M],[0 max(y)+max(y)/2],'color','green','LineWidth',2);
        line([lower(k) lower(k)],[0 max(y)-max(y)/3],'color','blue','LineWidth',2);
        line([upper(k) upper(k)],[0 max(y)-max(y)/3],'color','blue','LineWidth',2);
        title(['timepoint:',num2str(k)],...
            'FontWeight','Normal');
        set(gca,'ytick',[]);
        
        p(k+1,:) = get(h(k+1),'position');
    end
    
    
    height = p(3,2) + p(3,4) - p(k+1,2);
    width = p(k+1,1) + p(k+1,3) - p(3,1);
    h(k+2) = axes('position',[p(k,1) p(k,2) width height],'visible','off');
    h(k+2).XLabel.Visible = 'on';
    h(k+2).YLabel.Visible = 'on';
    
    axes(h(k+2));
    xlabel('Amplitude');
    ylabel({'Pobability Density';'at given IR time-point'});
    
    sgtitle({[mod_full,' vs Permutations (n=500)'];'LFP-LFP';['fmax = ',num2str(fmax(maxHz)),' Hz']});
    
    pos = get(gcf,'pos');
    set(gcf,'pos',[pos(1) pos(2) 10 10]);
    
    hold off
    
    % Paired t-test between wilson decomp. IR & mean of 500 permutations
    
    h12Nnew1=mean(h12N,2);
    t = ttest2(ir12s(1:30),h12Nnew1(1:30));
    if t==0
        fprintf('%s paired t-test for %u Hz is not significant\n',mod_full,fmax(maxHz));
    else
        fprintf('%s paired t-test for %u Hz is significant\n',mod_full,fmax(maxHz));
    end
end
