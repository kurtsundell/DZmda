











%Youngest Gaussian Fit of PDP (YGF)
x_YGF=0:0.1:4500;
pdp_YGF = pdp5(dist_data_perm(:,1),dist_data_perm(:,2),0,4500,0.1); %1 sigma pdp input data
pdp_YGF = pdp_YGF/1/sum(pdp_YGF); % normalize pdp to 1

[trs,trlocs] = findpeaks(-pdp_YGF,x_YGF);
if isempty(trlocs) == 0
	tridx = find(x_YGF==trlocs(1,1));	
else
	tridx = find(x_YGF==str2num(get(H.xmaxp,'string')));
end

[YGF_minr, YGF_minc] = find(pdp_YGF(1:tridx)>1E-6);
YGF_x = x_YGF(min(YGF_minc):tridx);
YGF_pdp = pdp_YGF(min(YGF_minc):tridx);

gofcells = cell(5,8);
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
AIC_values = zeros(8,2);
for k = 1:8
    if k == 1
        ft = fittype( 'gauss1' );
        [fitresult, gof] = fit( YGF_x', YGF_pdp', ft, opts );
        gofcells(:,k) = struct2cell(gof);
        Models(k,1:3*k) = coeffvalues(fitresult);  
        adj_r2(k) = 1 - ((1-cell2mat(gofcells(2,k)))*(length(YGF_x)-1)/(length(YGF_x)-3*k-1));      
    elseif k == 2
        ft = fittype( 'gauss2' ); 
        [fitresult, gof] = fit( YGF_x', YGF_pdp', ft, opts );
        gofcells(:,k) = struct2cell(gof);
        Models(k,1:3*k) = coeffvalues(fitresult);
        adj_r2(k) = 1 - ((1-cell2mat(gofcells(2,k)))*(length(YGF_x)-1)/(length(YGF_x)-3*k-1));
    elseif k == 3
        ft = fittype( 'gauss3' ); 
        [fitresult, gof] = fit( YGF_x', YGF_pdp', ft, opts );
        gofcells(:,k) = struct2cell(gof);
        Models(k,1:3*k) = coeffvalues(fitresult);
        adj_r2(k) = 1 - ((1-cell2mat(gofcells(2,k)))*(length(YGF_x)-1)/(length(YGF_x)-3*k-1));
    elseif k == 4
        ft = fittype( 'gauss4' );
        [fitresult, gof] = fit( YGF_x', YGF_pdp', ft, opts );
        gofcells(:,k) = struct2cell(gof);
        Models(k,1:3*k) = coeffvalues(fitresult);
        adj_r2(k) = 1 - ((1-cell2mat(gofcells(2,k)))*(length(YGF_x)-1)/(length(YGF_x)-3*k-1));
    elseif k == 5
        ft = fittype( 'gauss5' ); 
        [fitresult, gof] = fit( YGF_x', YGF_pdp', ft, opts );
        gofcells(:,k) = struct2cell(gof);
        Models(k,1:3*k) = coeffvalues(fitresult);
        adj_r2(k) = 1 - ((1-cell2mat(gofcells(2,k)))*(length(YGF_x)-1)/(length(YGF_x)-3*k-1));
    elseif k == 6
        ft = fittype( 'gauss6' ); 
        [fitresult, gof] = fit( YGF_x', YGF_pdp', ft, opts );
        gofcells(:,k) = struct2cell(gof);
        Models(k,1:3*k) = coeffvalues(fitresult);
        adj_r2(k) = 1 - ((1-cell2mat(gofcells(2,k)))*(length(YGF_x)-1)/(length(YGF_x)-3*k-1));
    elseif k == 7
        ft = fittype( 'gauss7' ); 
        [fitresult, gof] = fit( YGF_x', YGF_pdp', ft, opts );
        gofcells(:,k) = struct2cell(gof);
        Models(k,1:3*k) = coeffvalues(fitresult);
        adj_r2(k) = 1 - ((1-cell2mat(gofcells(2,k)))*(length(YGF_x)-1)/(length(YGF_x)-3*k-1));
    elseif k == 8
        ft = fittype( 'gauss8' ); 
        [fitresult, gof] = fit( YGF_x', YGF_pdp', ft, opts );
        gofcells(:,k) = struct2cell(gof);
        Models(k,1:3*k) = coeffvalues(fitresult);
        adj_r2(k) = 1 - ((1-cell2mat(gofcells(2,k)))*(length(YGF_x)-1)/(length(YGF_x)-3*k-1));
    end
end

[bestadjr2, bestfit] = find(diff(abs(cell2mat(gofcells(4,:))))<0);
bestfit=min(bestfit);
negcoef = 1;
while negcoef == 1
    coefficient = [];
    meanvals = [];
    sigvals = [];
    for k = 1:bestfit
    coefficient(k) = Models(bestfit, 3*k-2);
    meanvals(k) = Models(bestfit, 3*k-1);
    sigvals(k) = Models(bestfit, 3*k);
    end
    [ind1 ind2]=find(coefficient<=0); 
    if isempty(ind1) == 1
        negcoef = 0;
    else
        negcoef = 1;
        bestfit = bestfit-1;
    end
end            

[meanvals, I] = sort(meanvals);
coefficient = coefficient(I);
sigvals = sigvals(I);
for i=1:length(meanvals)
    if coefficient(i)>=0
        YGF = meanvals(i);
        YGF_1s = sigvals(i)/sqrt(2);
        break
    end
end

YGF_2s = YGF_1s*2;
x32 = [0:0.1:4400];
Yhat = normpdf(x32,YGF,YGF_1s);
Yhat = max(pdp)*(Yhat./max(Yhat));
% p1 = plot(x3,Yhat,'Color',[1 0 0]);
% set(p1,'linewidth',2)
% setheight = Yhat(round(YGF/0.1,0));
% s2 = scatter(YGF,setheight+spcy_s, 150, 'filled', 'v', 'markeredgecolor', 'k',  'markerfacecolor', [1 0 0], 'linewidth', 2);















%{


%%%%%%%%%%    Youngest Gaussian Fit Alternative       %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%truncate data at first negative inflection
slope = diff(pdp_YGF);
[slopeval,slopelocs1] = findpeaks(-slope);
[slopeind1 slopeind2] = find(slopeval>1E-10);

for i = 1:length(dist_data_perm(:,1))
if dist_data_perm(i,1) < x(slopelocs1(1,slopeind2)) 
YGF_data(i,:) = dist_data_perm(i,:);
else
YGF_data(i,1:2) = 0;
end
end
YGF_data = YGF_data(any(YGF_data ~= 0,2),:);
YGF_data = sortrows(YGF_data,1); 
%}


x_YGF=0:0.1:4500;
pdp_YGF = pdp5(dist_data_perm(:,1),dist_data_perm(:,2),0,4500,0.1); %1 sigma pdp input data
pdp_YGF = pdp_YGF/1/sum(pdp_YGF); % normalize pdp to 1

[trs,trlocs] = findpeaks(-pdp_YGF,x_YGF);


%truncate data at first PDP nadir
for i = 1:length(dist_data_perm(:,1))
if dist_data_perm(i,1) < trlocs(1,1)
YGF_data(i,:) = dist_data_perm(i,:);
else
YGF_data(i,1:2) = 0;
end
end
YGF_data = YGF_data(any(YGF_data ~= 0,2),:);
YGF_data = sortrows(YGF_data,1); 


%fit gaussians to ages comprising youngest PDP peak

for i=1:8
    if i>size(YGF_data,1)-1
        AICc = -99.99;
        break
    end
    try 
    GMM = fitgmdist(YGF_data(:,1),i,...
        'CovarianceType','full', 'SharedCov', false,...
        'Options',statset('MaxIter',10000),'ProbabilityTolerance',1E-10,'Replicates',10);
    AIC(i)=GMM.AIC;
    AICc(i) = AIC(i)+((2*2^2+2*2)/(length(YGF_data(:,1))-2-1));
    BIC(i)=GMM.BIC;
    catch
        break
    end
end
[min_peaks,min_int] = min(AICc);
if AICc == -99.99
    YGF = YGF_data(1,1);
    YGF_1s = YGF_data(1,2);
    YGF_2s = 2*YGF_1s;
else
    GMM = fitgmdist(YGF_data(:,1),min_int,...
        'CovarianceType','full','SharedCov', false,...
        'Options',statset('MaxIter',10000),'ProbabilityTolerance',1E-10,'Replicates',10);
    % for j=1:min_int
    %     line(x,pks(1,1)*normpdf(x,GMM.mu(j),sqrt(GMM.Sigma(j)))/...
    %         normpdf(GMM.mu(j),GMM.mu(j),sqrt(GMM.Sigma(j))),'color','green', 'linewidth',2);
    % end

    [YGF, min_int] = min(GMM.mu);
    YGF_1s = (GMM.Sigma(1,1,min_int))^0.5;
    YGF_2s = 2*YGF_1s;
end


x32 = [0:0.1:4400];
Yhat = normpdf(x32,YGF,YGF_1s);
Yhat = max(pdp)*(Yhat./max(Yhat));




%}

