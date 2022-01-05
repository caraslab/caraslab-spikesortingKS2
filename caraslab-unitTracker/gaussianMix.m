function gaussianMix(channel, unit, corrscore, wavescore, autoscore, basescore, survival)
% This function reclassifies the data to plot 
% the gaussian mixtures and decision boundary using a
% couple of metrics.
% Features plotted for illustration purposes are: 
% waveform similarity vs cross-correlation and 
% baseline firing vs waveform similarity.

sameChan = cell(size(corrscore));
sameUnit = cell(size(corrscore));
for iid=1:length(corrscore)
    sameChan{iid} = bsxfun(@eq,channel{iid}(:),channel{iid+1}(:)');
    sameUnit{iid} = sameChan{iid} & bsxfun(@eq,unit{iid}(:),unit{iid+1}(:)');
end

data = [unroll(corrscore) unroll(wavescore) unroll(autoscore) unroll(basescore)];
sameChan = logical(unroll(sameChan));
sameUnit = logical(unroll(sameUnit));
survival = logical(unroll(survival));

datamin = min(quantile(data(sameUnit,:),.01),quantile(data(~sameUnit,:),.01));
datamax = max(quantile(data(sameUnit,:),.99),quantile(data(~sameUnit,:),.99));
dataRange = [datamin; datamax] + [-.1; .1] * (datamax-datamin);
dataMean = mean([mean(data(survival,:)); mean(data(~survival,:))]);

negative = data(~sameChan,:);
[~,~,P] = classify(negative,data,survival,'quadratic');

negative = P(:,2);
% threshold = quantile(negative,.95);
threshold = quantile(negative,.95);
% threshold = .5;
fprintf('Threshold = %0.3f\n',threshold);

if(size(data,1)>2000)
    markerSizes = [.5 1.5];
    lineWidths = .2;
else
    markerSizes = [4 4];
    lineWidths = 1;
end

% Fit GMM
n = (2*threshold)/(1+2*threshold);
p = 1/(1+2*threshold);
prior = [n p];
[~,~,~,~,coeff] = classify(data,data,survival,'quadratic',prior);

meanSame = mean(data(survival,:));
meanDiff = mean(data(~survival,:));
covSame = cov(data(survival,:));
covDiff = cov(data(~survival,:));

contoursSame = mvnpdf(repmat(meanSame,4,1)+[sqrt(covSame(1,1)).*norminv(1-([.05 .25 .5 .75])./2)', zeros(4,3)], meanSame, covSame);
contoursDiff = mvnpdf(repmat(meanDiff,4,1)+[sqrt(covDiff(1,1)).*norminv(1-([.05 .25 .5 .75])./2)', zeros(4,3)], meanDiff, covDiff);

K = coeff(1,2).const;
L = coeff(1,2).linear; 
Q = coeff(1,2).quadratic;
Q(tril(true(4),-1)) = 0;

subplot(1,2,1);
hold on;
x_var_idx = 1;
y_var_idx = 2;
h = gscatter(data(:,x_var_idx),data(:,y_var_idx),survival,'br','.x',markerSizes,'off','Corr','Wave');
set(h,'LineWidth',lineWidths);
h2 = ezplot(@corrWaveBoundary, [-3 8 -3 8]);
set(h2,'Color','k','LineWidth',1,'LineStyle','-');

[X,Y] = meshgrid(range100(dataRange(:,x_var_idx)),range100(dataRange(:,y_var_idx)));
pSame = nan(size(X));
pSame(:) = mvnpdf([X(:),Y(:),repmat(meanSame(3:4),numel(X),1)], meanSame, covSame);
pSame = pSame ./ contoursSame(4);
contour(X,Y,pSame,contoursSame ./ contoursSame(4));
pDiff = nan(size(pSame));
pDiff(:) = mvnpdf([X(:),Y(:),repmat(meanDiff(3:4),numel(X),1)], meanDiff, covDiff);
pDiff = pDiff ./ contoursDiff(4);
contour(X,Y,pDiff,contoursDiff ./ contoursDiff(4));

xlim(dataRange(:,x_var_idx));
ylim(dataRange(:,y_var_idx));
title('');
xlabel('Pairwise Similarity');
ylabel('Wave Similarity');
set(gca,'XTick',atanh([0 .5 .9 .99 .999]));
set(gca,'XTickLabel',{'0','.5','.9','.99','.999'});
set(gca,'YTick',atanh([0 .5 .9 .99 .999]));
set(gca,'YTickLabel',{'0','.5','.9','.99','.999'});
axis square;


subplot(1,2,2);
x_var_idx = 4;
y_var_idx = 2;
hold on;
h = gscatter(data(:,x_var_idx),data(:,y_var_idx),survival,'br','.x',markerSizes,'off','Corr','Wave');
set(h,'LineWidth',lineWidths);
h2 = ezplot(@baseWaveBoundary, [-3 8 -3 8]);
set(h2,'Color','k','LineWidth',1,'LineStyle','-');

% MML edit: couldn't figure out why this code looks so different than previous
% plot 
[X,Y] = meshgrid(range100(dataRange(:,x_var_idx)),range100(dataRange(:,y_var_idx)));
pSame = nan(size(X));
pSame(:) = mvnpdf([repmat(meanSame(1),numel(X),1),Y(:),repmat(meanSame(3),numel(X),1),X(:)], meanSame, covSame);
pSame = pSame ./ contoursSame(4);
contour(X,Y,pSame,contoursSame ./ contoursSame(4));
pDiff = nan(size(pSame));
pDiff(:) = mvnpdf([repmat(meanDiff(1),numel(X),1),Y(:),repmat(meanDiff(3),numel(X),1),X(:)], meanDiff, covDiff);
pDiff = pDiff ./ contoursDiff(4);
contour(X,Y,pDiff,contoursDiff ./ contoursDiff(4));

xlim(dataRange(:,x_var_idx));
ylim(dataRange(:,y_var_idx));
title('');
xlabel('\Delta Log Mean Rate');
ylabel('Wave Similarity');
set(gca,'XTick',-3:3);
set(gca,'XTickLabel',-3:3);
set(gca,'YTick',atanh([0 .5 .9 .99 .999]));
set(gca,'YTickLabel',{'0','.5','.9','.99','.999'});
axis square;

% subplot(2,2,3);
% hold on;
% h = gscatter(data(:,4),data(:,2),survival,'br','.x',markerSizes,'off','Corr','Wave');
% set(h,'LineWidth',lineWidths);
% h2 = ezplot(@baseWaveBoundary, [-3 8 -3 8]);
% set(h2,'Color','k','LineWidth',1,'LineStyle','-');
% 
% [X,Y] = meshgrid(range100(dataRange(:,4)),range100(dataRange(:,2)));
% pSame = nan(size(X));
% pSame(:) = mvnpdf([repmat(meanSame(1),numel(X),1),Y(:),repmat(meanSame(3),numel(X),1),X(:)], meanSame, covSame);
% pSame = pSame ./ contoursSame(4);
% contour(X,Y,pSame,contoursSame ./ contoursSame(4));
% pDiff = nan(size(pSame));
% pDiff(:) = mvnpdf([repmat(meanDiff(1),numel(X),1),Y(:),repmat(meanDiff(3),numel(X),1),X(:)], meanDiff, covDiff);
% pDiff = pDiff ./ contoursDiff(4);
% contour(X,Y,pDiff,contoursDiff ./ contoursDiff(4));
% 
% xlim(dataRange(:,4));
% ylim(dataRange(:,2));
% title('');
% xlabel('\Delta Log Mean Rate');
% ylabel('Wave Similarity');
% set(gca,'XTick',-3:3);
% set(gca,'XTickLabel',-3:3);
% set(gca,'YTick',atanh([0 .5 .9 .99 .999]));
% set(gca,'YTickLabel',{'0','.5','.9','.99','.999'});
% axis square;


% keyboard;
% 
% % XXXXX This uses gaussian data to assess accuracy.  Needs to switch to
% % integrating only over data points
testData = [mvnrnd(meanSame,covSame,5000); mvnrnd(meanDiff,covDiff,5000)];
% testClass = [ones(5000,1); repmat(2,5000,1)];
% 
% sensitivityTable = nan(4);
% specificityTable = nan(4);
% for i=1:4
%     cp = testVars(i);
%     specificityTable(i,i) = cp.Specificity;
%     sensitivityTable(i,i) = cp.Sensitivity;
% end
% combos = nchoosek(1:4,2);
% for i=1:size(combos,1)
%     cp = testVars(combos(i,:));
%     specificityTable(combos(i,1),combos(i,2)) = cp.Specificity;
%     sensitivityTable(combos(i,1),combos(i,2)) = cp.Sensitivity;
% end
% cp = testVars(1:4);
% 
% probabilityTrue = mean(C(sameChan));
% correct = probabilityTrue*cp.Sensitivity + (1-probabilityTrue)*cp.Specificity;

set(gcf,'PaperPosition', [0 0 4.5 4.5]);

function cp = testVars(vars)
classOut = classify(testData(:,vars), data(:,vars), survival, 'quadratic');
classOut = double(classOut);
classOut(classOut==0) = 2;
cp = classperf(testClass, classOut);
end

function y = corrWaveBoundary(corr, wave)
x = dataMean;
x(1) = corr;
x(2) = wave;

y = K + x*L + sum(sum((x'*x).*Q));
end

function y = baseWaveBoundary(base, wave)
x = dataMean;
x(4) = base;
x(2) = wave;

y = K + x*L + sum(sum((x'*x).*Q));
end

function y = range100(r)
    y = r(1):range(r)/100:r(2);
end

end