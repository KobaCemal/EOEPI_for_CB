clear

ls=dir('cleanrs/sc');
list={ls.name};
list=list(3:end);
path=ls.folder;

matrices=nan(132,400,28);


for i=1:14
    matrices(:,:,i)=readmatrix([path '/' list{i}]);
end


ls=dir('wav/sc');
list={ls.name};
list=list(3:end);
path=ls.folder;
for i=1:14
    matrices(:,:,i+14)=readmatrix([path '/' list{i}]);
end

%% Comparison of the network metrics
matrix_of_interest=matrices;

numsub=28;
numrois=400;
correlation_matrices=nan(numrois,numrois,28);
transformed_matrices=nan(numrois,numrois,28);
for i=1:28
    correlation_matrices(:,:,i)=corr(matrix_of_interest(:,:,i));
    transformed_matrices(:,:,i)=atanh(correlation_matrices(:,:,i));
end

spar=[0.1 0.2 0.3];
statsAll01=nan(numsub,18);
statsAll02=nan(numsub,18);
statsAll03=nan(numsub,18);

for i=1:numsub
    disp(i)
    cor=correlation_matrices(:,:,i);%% select the correlation matrix
    parfor (j=1:length(spar),6)
        thr=threshold_proportional(cor,spar(j));%% threshold it with the sparsity level
        bin=weight_conversion(thr, 'binarize');
        len=weight_conversion(thr, 'lengths');
        nor=weight_conversion(thr, 'normalize');
        
        deg=degrees_und(thr);%%node degree
        %     for j =1:500
        %         if deg(j)>mean(deg)+std(deg);
        %             hubwith(j,1)=1;
        %         else
        %             hubwith(j,1)=0;
        %         end
        %     end

        str=strengths_und(thr);%%strengths
        %     for j =1:500
        %         if str(j)>mean(str)+std(str);
        %             hubwith(j,2)=1;
        %         else
        %             hubwith(j,2)=0;
        %         end
        %     end

        den=density_und(thr);%%density

        clu=clustering_coef_wu(nor);%%cluster coef

        tra=transitivity_wu(nor);%%transitivity

        ass=assortativity_wei(thr, 0);%%assortativity

        eff=efficiency_wei(thr);%%global efficiency

        [Ci,Q] = modularity_und(thr);

        bet=betweenness_wei(len);%betweenness centrality

        sparsitylevel(j,:)=[max(deg) min(deg(deg>0)) mean(deg(deg>0)) ...
            max(str) min(str(str>0)) mean(str(str>0)) ...
            den ...
            max(clu) min(clu(clu>0)) mean(str(clu>0)) ...
            tra ...
            ass ...
            eff ...
            max(Ci) Q ...
            max(bet) min(bet(bet>0)) mean(bet(bet>0))];


        %     for j =1:500
        %         if bet(j)>mean(bet)+std(bet);
        %             hubwith(j,3)=1;
        %         else
        %             hubwith(j,3)=0;
        %         end
        %     end

        % eff2=efficiency_wei(thr,2);%%local efficiency
        %     for i =1:500
        %         if deg(i)>mean(eff2)+std(eff2);
        %             hubwith(i,4)=1;
        %         else
        %             hubwith(i,4)=0;
        %         end
        %     end

        % for j=1:500;  hubwith(j,4)=hubwith(j,1)+hubwith(j,2)+hubwith(j,3);end


        % hubwithAll(:,i)=hubwith(:,4);
    end
    statsAll01(i,:)=sparsitylevel(1,:);
    statsAll02(i,:)=sparsitylevel(2,:);
    statsAll03(i,:)=sparsitylevel(3,:);
end


for i=1:min(size(statsAll01))
    [h,p,ci,stats] = ttest(statsAll01(1:14,i),statsAll01(15:28,i));
    tp01(1,i)=(mean(statsAll01(1:14,i))-mean(statsAll01(15:28,i)))/abs(mean(statsAll01(1:14,i)))*100;
    tp01(2,i)=(mean(statsAll01(1:14,i))-mean(statsAll01(15:28,i)));
    tp01(3,i)=h;
    tp01(4,i)=p;
    tp01(5,i)=stats.tstat;
    tp01(6,i)=abs(computeCohen_d(statsAll01(1:14,i),statsAll01(15:28,i), 'paired'));

end
tp01=tp01'

for i =1:min(size(statsAll02))
    [h,p,ci,stats] = ttest(statsAll02(1:14,i),statsAll02(15:28,i));
    tp02(1,i)=(mean(statsAll02(1:14,i))-mean(statsAll02(15:28,i)))/abs(mean(statsAll02(1:14,i)))*100;
    tp02(2,i)=(mean(statsAll02(1:14,i))-mean(statsAll02(15:28,i)));
    tp02(3,i)=h;
    tp02(4,i)=p;
    tp02(5,i)=stats.tstat;
    tp02(6,i)=abs(computeCohen_d(statsAll02(1:14,i),statsAll02(15:28,i), 'paired'));
end
tp02=tp02'

for i =1:min(size(statsAll03))
    [h,p,ci,stats] = ttest(statsAll03(1:14,i),statsAll03(15:28,i));
    tp03(1,i)=(mean(statsAll03(1:14,i))-mean(statsAll03(15:28,i)))/abs(mean(statsAll03(1:14,i)))*100;
    tp03(2,i)=(mean(statsAll03(1:14,i))-mean(statsAll03(15:28,i)));
    tp03(3,i)=h;
    tp03(4,i)=p;
    tp03(5,i)=stats.tstat;
    tp03(6,i)=abs(computeCohen_d(statsAll03(1:14,i),statsAll03(15:28,i), 'paired'));
end
tp03=tp03'

%% Node-wise comparison of the matrices
clear

ls=dir('cleanrs/sc');
list={ls.name};
list=list(3:end);
path=ls.folder;

matrices=nan(132,400,28);


for i=1:14
    matrices(:,:,i)=readmatrix([path '/' list{i}]);
end


ls=dir('wav/sc');
list={ls.name};
list=list(3:end);
path=ls.folder;
for i=1:14
    matrices(:,:,i+14)=readmatrix([path '/' list{i}]);
end

matrix_of_interest=matrices;

numsub=28;
numrois=400;
correlation_matrices=nan(numrois,numrois,28);
transformed_matrices=nan(numrois,numrois,28);
for i=1:28
    correlation_matrices(:,:,i)=corr(matrix_of_interest(:,:,i));
    transformed_matrices(:,:,i)=atanh(correlation_matrices(:,:,i));
end

ps=nan(numrois,numrois);
ts=nan(numrois,numrois);
ds=nan(numrois,numrois);
for i=1:numrois
    disp(i)
    for j=1:numrois
        [h,p,ci,stats] = ttest(squeeze(transformed_matrices(i,j,1:14)),squeeze(transformed_matrices(i,j,15:28)));
        ps(i,j)=p;
        ts(i,j)=stats.tstat;
%         ds(i,j)=computeCohen_d(squeeze(transformed_matrices(i,j,1:14)),squeeze(transformed_matrices(i,j,15:28)),'paired');
    end
end
vectorized_ps=reshape(triu(ps),[],1);
fdr_bh(vectorized_ps(vectorized_ps>0));
imagesc(ts.*(abs(ts)>3))
% on thresholded matrices
% for i=1:28
%     cor=correlation_matrices(:,:,i);
%     thr(:,:,i)=threshold_proportional(cor,0.05);
% end
% 
% ps=nan(numrois,numrois);
% ts=nan(numrois,numrois);
% ds=nan(numrois,numrois);
% for i=1:numrois
%     disp(i)
%     for j=1:numrois
%         [h,p,ci,stats] = ttest(squeeze(thr(i,j,1:14)),squeeze(thr(i,j,15:28)));
%         ps(i,j)=p;
%         ts(i,j)=stats.tstat;
%         ds(i,j)=computeCohen_d(squeeze(thr(i,j,1:14)),squeeze(thr(i,j,15:28)),'paired');
%     end
% end
% vectorized_ps=reshape(triu(ps),[],1);
% fdr_bh(vectorized_ps(vectorized_ps>0),0.05,'dep','yes');

%% Difference between visual and somatomotor

vis=squeeze(matrices(:,221,:));
som=squeeze(matrices(:,265,:));

for i=1:14
    corr_before(i)=atanh(corr(vis(:,i),som(:,i), type='Spearman'));
    corr_after(i)=atanh(corr(vis(:,i+14),som(:,i+14), type='Spearman'));
    %corr_eoepi_vis(i)=corr(vis(:,i),eoepi(:,i));
    %corr_eoepi_som(i)=corr(som(:,i),eoepi(:,i));
end

mean(corr_before)
mean(corr_after)

std(corr_before)
std(corr_after)

[H2,P2,CI2,stats2] = ttest(corr_before,corr_after)

vars=var(eoepi)
corr(vars',corr_before',type='Spearman')
corr(vars',corr_after', type='Spearman')




