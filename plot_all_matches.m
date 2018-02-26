%plot quantification results for matched data across different samples
%input matrix as follows:
% rows: sorted Pr (release probability) quartiles for each dataset
% columns: Pr Pr(std) Brp Brp(std) Cac Cac(std)

color_mat = [0.5 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 1 0.5 0; 1 0 0.5;];
plot_num = length(wt_matches)/4;


cat_select = 7; %7 for dist
cat_select_2 = 5;

mean_counts_total = mean(wt_matches(:,cat_select));

mean_pr_total = mean(wt_matches(:,cat_select_2));


MarkerSize = 5;
   clf
   for i=1:plot_num
        hold on
%         for q=1:num_bins
        idx_start = 4*(i-1)+1;
        idx_end = 4*(i);
        
        data_current = wt_matches(idx_start:idx_end,:);
        s = mean(data_current(:,cat_select))-mean_counts_total;
        s1 = mean(data_current(:,cat_select_2))-mean_pr_total;
        brp_norm(i,:) = data_current(:,cat_select)-s;
        pr_binned(i,:) = data_current(:,cat_select_2);
        pr_norm(i,:) = data_current(:,cat_select_2)-s1;
        pr_std(i,:) = data_current(:,cat_select_2+1);
        brp_std(i,:) = data_current(:,cat_select+1);
        
       
   end
   %%
   for i=1:plot_num
       
        hold on
        plot(brp_norm(i,:),pr_norm(i,:),'Color',color_mat(i,:),'Marker','o')
        errorbarxy(brp_norm(i,:),pr_norm(i,:),brp_std(i,:),pr_std(i,:),color_mat(i,:));

        
   end 
   
       %for mean calculations
        brp_mean = mean(brp_norm,1);
        pr_mean = mean(pr_norm,1);
        [b1,b2,r1,r2] = linreg(brp_mean',pr_mean');
      
%         clf
%         
        hold on
%         xline = [1 0; 1 1];
%         yline = xline*b2;
%         line([0;1],yline,'LineWidth',2)
        plot(brp_mean,pr_mean,'k.','MarkerSize',30)
% %% fitting
% modelFun = @(b,x) b(1).*(1-exp(-b(2).*x))
% start = [650; 25]
% wnlm = fitnlm(pr_mean,brp_mean,modelFun,start)
% xx = linspace(0,0.5)';
% line(xx,predict(wnlm,xx),'linestyle','--','color','k')
%       set(leg,'location','NorthWest')
        XLimits=xlim;
        YLimits=ylim;
        xlim([0,max(brp_mean)*1.2])
        ylim([0,max(pr_mean)*1.8])
        xlabel('Distance to edge (px)')
        ylabel('Cac')
        FigureStandardizer(get(gca,'xlabel'), get(gca,'ylabel'), gca);
        set(gcf, 'color', 'white');
%         set(leg,'location','EastOutside')
        set(gcf,'units','normalized','position',[0.1,0.1,0.6,0.8])
        set(gca,'position',[0.15,0.25,0.5,0.5])