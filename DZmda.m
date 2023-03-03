% DZMDA MATLAB code for DZmda.fig %%

% SET DEFAULT COMMAND LINE AND HANDLE STRUCTURE %%
function varargout = DZmda(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',mfilename,'gui_Singleton',gui_Singleton,'gui_OpeningFcn',@DZmda_OpeningFcn,'gui_OutputFcn',@DZmda_OutputFcn,'gui_LayoutFcn',[],'gui_Callback',[]);
if nargin && ischar(varargin{1})
	gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
	[varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
	gui_mainfcn(gui_State, varargin{:});
end
function DZmda_OpeningFcn(hObject, eventdata, H, varargin)
H.output = hObject;
guidata(hObject, H);
function varargout = DZmda_OutputFcn(hObject, eventdata, H)
varargout{1} = H.output;
set(H.input1s,'Value',1)
H.export_dist = 0;
%plot_distribution(hObject, eventdata, H)
guidata(hObject,H);

% DISTRIBUTION PLOTTER %%
function load_Callback(hObject, eventdata, H)
[filename pathname] = uigetfile({'*'},'File Selector');
data = readtable(char(strcat(pathname, filename)));
data = table2array(data);
data = sortrows(data);
set(H.uitable1, 'Data', data);
set(H.setx,'Value',1)
%data = get(H.uitable1, 'Data');
%set(H.filepath,'String',[strcat(pathname, filename)])
plot_distribution(hObject, eventdata, H)

function plot_distribution(hObject, eventdata, H)

if H.export_dist == 1
	figure;
end
if H.export_dist == 0
	cla(H.axes_comp,'reset');
	axes(H.axes_comp);
end
H.export_dist = 0;
guidata(hObject,H);

data = get(H.uitable1, 'Data');

if iscell(data) == 1
	data = cell2num(data);
end

if get(H.input1s,'Value') == 1
	columnname = {'Age', '±1s'};
	set(H.uitable1,'ColumnName', columnname);
	dist_data = data;
end

if get(H.input2s,'Value') == 1
	columnname = {'Age', '±2s'};
	set(H.uitable1, 'ColumnName', columnname);
	dist_data(:,1) = data(:,1);
	dist_data(:,2) = data(:,2)./2;
end

n = length(dist_data(:,1));

%%%% All calculations assume 1 sigma % input below. Adjusted above if needed.

hold on

%keep original data
dist_data_perm = dist_data;

% for i = 1:length(dist_data(:,1))
% 	if dist_data(i,1) > str2double(get(H.xminp,'String')) && dist_data(i,1) < str2double(get(H.xmaxp,'String'))
% 		dist_data(i,:) = dist_data(i,:);
% 	else
% 		dist_data(i,1:2) = 0;
% 	end
% end
% 
% dist_data = dist_data(any(dist_data ~= 0,2),:);
dist_data = sortrows(dist_data,1);

xint = 0.1;

x=0:xint:4500;

% Make PDP
pdp = pdp5(dist_data(:,1),dist_data(:,2),0,4500,xint); %1 sigma pdp input data
pdp = pdp/1/sum(pdp); % normalize pdp to 1


% MDA calculations
spcy_s = 0.0003; %adjust scatter triangles up
%spcx_t = 0; %adjust text to left
%spcy_t = 0.001; %adjust text up

% Youngest single grain (YSG)
YSG_hi_tmp = dist_data(:,1) + dist_data(:,2);
[~,sortIdx0] = sort(YSG_hi_tmp);
dist_data_sort0 = dist_data(sortIdx0,:);
YSG = dist_data_sort0(1,1);
YSG_1s = dist_data_sort0(1,2);
YSG_2s = 2*dist_data_sort0(1,2);

% Youngest graphical peak of PDP (YPP)
[pks,locs] = findpeaks(pdp,x);

%text(locs(1,1)-spcx_t,pks(1,1)+spcy_t,num2str(locs(1,1)), 'FontSize',16, 'Color', 'b', 'horizontalAlignment', 'left')
YPP = locs(1,1);

























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
f = fit(YGF_x',YGF_pdp','gauss1');
YGF = f.b1;
YGF_1s = f.c1/sqrt(2);
YGF_2s = YGF_1s*2;
x32 = [0:0.1:4500];
Yhat = feval(f,x32);









%{



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






%}







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
    YGF3 = YGF_data(1,1);
    YGF3_1s = YGF_data(1,2);
    YGF3_2s = 2*YGF3_1s;
else
    GMM = fitgmdist(YGF_data(:,1),min_int,...
        'CovarianceType','full','SharedCov', false,...
        'Options',statset('MaxIter',10000),'ProbabilityTolerance',1E-10,'Replicates',10);
    % for j=1:min_int
    %     line(x,pks(1,1)*normpdf(x,GMM.mu(j),sqrt(GMM.Sigma(j)))/...
    %         normpdf(GMM.mu(j),GMM.mu(j),sqrt(GMM.Sigma(j))),'color','green', 'linewidth',2);
    % end

    [YGF3, min_int] = min(GMM.mu);
    YGF3_1s = (GMM.Sigma(1,1,min_int))^0.5;
    YGF3_2s = 2*YGF3_1s;
end



%}















if get(H.setx,'Value') == 1
	
	xmin = round(YGF - 10*YGF_1s);
	xmax = round(YGF + 10*YGF_1s);
	
% 	xmin = round(min(dist_data(:,1)-2*dist_data(:,2)) - min(dist_data(:,1)-dist_data(:,2))*.05); % make nice plots
% 	xmax = round(max(dist_data(:,1)+2*dist_data(:,2)) +  max(dist_data(:,1)+dist_data(:,2)).*.05); % make nice plots
	set(H.xminp,'String',xmin)
	set(H.xmaxp,'String',xmax)
end

if get(H.setx,'Value') == 0
	xmin = str2num(get(H.xminp,'String'));
	xmax = str2num(get(H.xmaxp,'String'));
end













% Youngest grain cluster at 1s (YGC1s)
YGC1s_hi_tmp = dist_data(:,1) + dist_data(:,2);
[~,sortIdx1] = sort(YGC1s_hi_tmp);
dist_data_sort1 = dist_data(sortIdx1,:);
YGC1s_hi = dist_data_sort1(:,1) + dist_data_sort1(:,2);
YGC1s_lo = dist_data_sort1(:,1) - dist_data_sort1(:,2);
for i = 1:length(dist_data_sort1(:,1))
	YGC1s_lo_LT(i,1) = YGC1s_lo(i,1) < YGC1s_hi(1,1);
end
for i = 1:length(YGC1s_lo_LT)
	if YGC1s_lo_LT(i,1) == 0
		break
	end
	YGC1s_data(i,:) = dist_data_sort1(i,:);
end
[YGC1s,YGC1s_1s,YGC1s_2s,YGC1s_mswd] = wm(YGC1s_data);

% Youngest grain cluster at 2s (YGC2s)
YGC2s_hi_tmp = dist_data(:,1) + 2*dist_data(:,2);
[~,sortIdx2] = sort(YGC2s_hi_tmp);
dist_data_sort2 = dist_data(sortIdx2,:);
YGC2s_hi = dist_data_sort2(:,1) + 2*dist_data_sort2(:,2);
YGC2s_lo = dist_data_sort2(:,1) - 2*dist_data_sort2(:,2);
for i = 1:length(dist_data_sort2(:,1))
	YGC2s_lo_LT(i,1) = YGC2s_lo(i,1) < YGC2s_hi(1,1);
end
for i = 1:length(YGC2s_lo_LT)
	if YGC2s_lo_LT(i,1) == 0
		break
	end
	YGC2s_data(i,:) = dist_data_sort2(i,:);
end
[YGC2s,YGC2s_1s,YGC2s_2s,YGC2s_mswd] = wm(YGC2s_data);

% Youngest three zircons (Y3Zo and Y3Za)
Y3Zo_sort = sortrows(YGC2s_data);
if length(Y3Zo_sort(:,1)) >= 3
	Y3Zo_data = Y3Zo_sort(1:3,:);
	[Y3Zo,Y3Zo_1s,Y3Zo_2s,Y3Zo_mswd]  = wm(Y3Zo_data);
else
	Y3Zo_data = Y3Zo_sort;
	Y3Zo = 0;
	Y3Zo_1s = 0;
	Y3Zo_2s = 0;
	Y3Zo_mswd = 0;
end
Y3Za_data = dist_data(1:3,:);
[Y3Za,Y3Za_1s,Y3Za_2s,Y3Za_mswd]  = wm(Y3Za_data);


% Tau method
[~,locsT] = findpeaks(-pdp);
if isempty(locsT) == 1
	TAU_lim = xmax;
end
if isempty(locsT) == 0
	TAU_lim = x(locsT(1,1));
end
TAU_data = dist_data(dist_data(:,1)<TAU_lim,:);
if length(TAU_data(:,1)) >= 3
	[TAU,TAU_1s,TAU_2s,TAU_mswd]  = wm(TAU_data);
end

% Youngest Statistical Population (YSP)
for i = 1:length(dist_data(:,1))-1
	YSP_data = dist_data(1:i+1,:);
	[YSP_wm(i,1),YSP_1s(i,1),YSP_2s(i,1),YSP_mswd(i,1)] = wm(YSP_data);
end
YSP_sub = 1 - YSP_mswd;
[~,YSP_idx] = min(abs((YSP_sub)));
YSP = YSP_wm(YSP_idx,1);
YSP_1s = YSP_1s(YSP_idx,1);
YSP_2s = YSP_2s(YSP_idx,1);
YSP_mswd = YSP_mswd(YSP_idx,1);


% ylabel('Frequency')
% %set(gca,'YTickLabel',[])
% set(gca,'XTickLabel',[])
% set(gca, 'box', 'on')
% legend('YDZ Monte Carlo Model','fontsize',14)
% legend boxoff
% set(gca, 'fontsize',12)













axes(H.axes_comp);

%gauss_adj = Yhat*(1/(max(Yhat)/pks(1,1)));

%gauss_adj1 = Yhat1*(1/(max(Yhat1)/pks(1,1)));
gauss_adj = Yhat*(1/(max(Yhat)/pks(1,1)));

p = plot(x, pdp, 'Color', 'b', 'LineWidth', 2);
p1 = plot(x32,gauss_adj,'Color',[1 0 0]);
%p2 = plot(x32,gauss_adj,'Color','g');

% for j=1:1
%     p3 = line(x,pks(1,1)*normpdf(x,GMM.mu(j),sqrt(GMM.Sigma(j)))/normpdf(GMM.mu(j),GMM.mu(j),sqrt(GMM.Sigma(j))),'color','m', 'linewidth',2);
% end


set(p,'linewidth',2)
set(p1,'linewidth',2)
%set(p2,'linewidth',2)


s1 = scatter(locs(1,1),pks(1,1), 150, 'filled', 'v', 'markeredgecolor', 'k',  'markerfacecolor', 'b', 'linewidth', 2);
s2 = scatter(YGF,pks(1,1), 150, 'filled', 'v', 'markeredgecolor', 'k',  'markerfacecolor', [1 0 0], 'linewidth', 2);





lgnd = legend([p, p1, s1, s2], 'Probability Density Plot', 'Youngest Gaussian Fit', 'YPP', 'YGF', 'fontsize',12);
set(lgnd,'Color','w');
set(gca, 'box', 'on')
%legend boxoff
set(gca, 'fontsize',14)

xlabel('Age (Ma)','Color','k')
ylabel('Probability','Color','k')
xlim([xmin xmax])














% Make Age plot
cla(H.axes_wm,'reset');
axes(H.axes_wm);

hold on % hold the line, duh nuh nuh nuh nuh..... love isn't always on time





if get(H.cluster_sort,'Value') == 0
	
	data2 = dist_data(dist_data(:,1) < xmax,:);
	data2(:,2) = data2(:,2).*2;
	len = length(data2(:,1));
	
end

if get(H.cluster_sort,'Value') == 1
	
	data2_hi_tmp = dist_data(:,1) + dist_data(:,2);
	[~,sortIdx2] = sort(data2_hi_tmp);
	data2 = dist_data(sortIdx2,:);
	
	
	data2 = data2(data2(:,1) < xmax,:);
	len = length(data2(:,1));
end









x2 = 1:1:len;

xmin2 = 0; % make nice plots
xmax2 = len+1; % make nice plots


u1 = plot([x2; x2], [(data2(:,1)+data2(:,2))'; (data2(:,1)-data2(:,2))'], '-r', 'Color', [.4 .6 1], 'LineWidth',5); % Error bars, much nicer than the errorbar function
u2 = plot([x2; x2], [(data2(:,1)+data2(:,2))'; (data2(:,1)-data2(:,2))'], '-r', 'Color', 'b', 'LineWidth',5); % Error bars, much nicer than the errorbar function

scatter(x2, data2(:,1), 175, 'k', 'filled','d')

hold off














xlim([xmin2 xmax2])
ylim([xmin xmax])

set(gca,'YTickLabel',[])
%set(gca,'XTickLabel',[])
xlabel('Rank Number')
set(gca, 'box', 'on')
set(gca, 'fontsize',14)

%ffff = legend([u1, u2], '± 2σ', '± 1σ');



camroll(-90) % rotate 90 deg
















% Make results plot
cla(H.axes_mda,'reset');
axes(H.axes_mda);

hold on

% Compile results
if length(YGC1s_data(:,1)) < 2
	YGC1s = 0;
	YGC1s_1s = 0;
	YGC1s_2s = 0;
	YGC1s_mswd = 0;
end
if length(YGC2s_data(:,1)) < 3
	YGC2s = 0;
	YGC2s_1s = 0;
	YGC2s_2s = 0;
	YGC2s_mswd = 0;
end
if length(TAU_data(:,1)) < 3
	TAU = 0;
	TAU_1s = 0;
	TAU_2s = 0;
	TAU_mswd = 0;
end

X = [1:1:9];
Y = [YSG;YPP;YGF;YGC1s;YGC2s;Y3Zo;Y3Za;TAU;YSP];
Y_hi = [0.0001;0.0001;0.0001;YGC1s_2s;YGC2s_2s;Y3Zo_2s;Y3Za_2s;TAU_2s;YSP_2s];
Y_lo = [0.0001;0.0001;0.0001;YGC1s_2s;YGC2s_2s;Y3Zo_2s;Y3Za_2s;TAU_2s;YSP_2s];

if get(H.sety,'Value') == 1
	yminp = round(min(nonzeros(Y(:,1))-nonzeros(Y_lo(:,1))) - 1); % make nice plots
	ymaxp = round(max(nonzeros(Y(:,1))+nonzeros(Y_hi(:,1))) + 1); % make nice plots
	set(H.ymin,'String',yminp)
	set(H.ymax,'String',ymaxp)
end

if get(H.sety,'Value') == 0
	yminp = str2num(get(H.ymin,'String'));
	ymaxp = str2num(get(H.ymax,'String'));
end



x3 = [1:1:9];
plot([x3; x3], [(Y+Y_hi)'; (Y-Y_lo)'], '-r', 'Color', [.4 .6 1], 'LineWidth',7) % Error bars, much nicer than the errorbar function
plot([x3; x3], [(Y+Y_hi./2)'; (Y-Y_lo./2)'], '-r', 'Color', 'k', 'LineWidth',7) % Error bars, much nicer than the errorbar function
scatter(X, Y, 250, 'filled','s', 'markeredgecolor','k','markerfacecolor','w','linewidth',2)
ylabel('Age (Ma)')
xlim([0 11])
ylim([yminp ymaxp])
set(gca, 'box', 'on')
names = {'YSG'; 'YPP'; 'YGF'; 'YGC1s'; 'YGC2s'; 'Y3Zo'; 'Y3Za'; 'YDZ'; 'TAU'; 'YSP'};
set(gca,'xtick',[1:10],'xticklabel',names)
set(gca, 'fontsize',14)

xlabel('Method')

set(gca,'YDir','reverse')




























tdata = [{sprintf('%1.2f',YSG)},{sprintf('%1.2f',YSG_1s)},{sprintf('%1.2f',YSG_2s)},{sprintf('%.f',1)},{'NA'};
	{sprintf('%1.2f',YPP)},{'NA'},{'NA'},{'NA'},{'NA'};
	{sprintf('%1.2f',YGF)},{sprintf('%1.2f',YGF_1s)},{sprintf('%1.2f',YGF_2s)},{'NA'},{'NA'};
	{sprintf('%1.2f',YGC1s)},{sprintf('%1.2f',YGC1s_1s)},{sprintf('%1.2f',YGC1s_2s)},{num2str(length(YGC1s_data(:,1)))},{sprintf('%1.2f',YGC1s_mswd)};
	{sprintf('%1.2f',YGC2s)},{sprintf('%1.2f',YGC2s_1s)},{sprintf('%1.2f',YGC2s_2s)},{num2str(length(YGC2s_data(:,1)))},{sprintf('%1.2f',YGC2s_mswd)};
	{sprintf('%1.2f',Y3Zo)},{sprintf('%1.2f',Y3Zo_1s)},{sprintf('%1.2f',Y3Zo_2s)},{num2str(length(Y3Zo_data(:,1)))},{sprintf('%1.2f',Y3Zo_mswd)};
	{sprintf('%1.2f',Y3Za)},{sprintf('%1.2f',Y3Za_1s)},{sprintf('%1.2f',Y3Za_2s)},{num2str(length(Y3Za_data(:,1)))},{sprintf('%1.2f',Y3Za_mswd)};
	{sprintf('%1.2f',TAU)},{sprintf('%1.2f',TAU_1s)},{sprintf('%1.2f',TAU_2s)},{num2str(length(TAU_data(:,1)))},{sprintf('%1.2f',TAU_mswd)};
	{sprintf('%1.2f',YSP)},{sprintf('%1.2f',YSP_1s)},{sprintf('%1.2f',YSP_2s)},{num2str(YSP_idx+1)},{sprintf('%1.2f',YSP_mswd)}];

if length(YGC1s_data(:,1)) < 2
	tdata(4,1:3) = {'n < 2'};
	tdata(4,5) = {'n < 2'};
end

if length(YGC2s_data(:,1)) < 3
	tdata(5,1:3) = {'n < 3'};
	tdata(5,5) = {'n < 3'};
end

if Y3Zo == 0
	tdata(6,1:3) = {'n < 3'};
	tdata(6,5) = {'n < 3'};
end

if length(TAU_data(:,1)) < 3
	tdata(8,1:3) = {'n < 3'};
	tdata(8,5) = {'n < 3'};
end

set(H.uitable6,'data',tdata)





function xmin_Callback(hObject, eventdata, H)
plot_distribution(hObject, eventdata, H)

function xmax_Callback(hObject, eventdata, H)
plot_distribution(hObject, eventdata, H)

function xint_Callback(hObject, eventdata, H)
plot_distribution(hObject, eventdata, H)

function uitable1_CellEditCallback(hObject, eventdata, H)
plot_distribution(hObject, eventdata, H)

function pastetable_Callback(hObject, eventdata, H)
data = paste;
set(H.uitable1, 'Data', data);
plot_distribution(hObject, eventdata, H)

function input1s_Callback(hObject, eventdata, H)
set(H.input1s,'Value',1)
set(H.input2s,'Value',0)
plot_distribution(hObject, eventdata, H)

function input2s_Callback(hObject, eventdata, H)
set(H.input1s,'Value',0)
set(H.input2s,'Value',1)
plot_distribution(hObject, eventdata, H)

function setx_Callback(hObject, eventdata, H)
plot_distribution(hObject, eventdata, H)

function xminp_Callback(hObject, eventdata, H)
set(H.setx,'Value',0)
plot_distribution(hObject, eventdata, H)

function xmaxp_Callback(hObject, eventdata, H)
set(H.setx,'Value',0)
plot_distribution(hObject, eventdata, H)

function exportplot_Callback(hObject, eventdata, H)
H.export_dist1 = 1;
H.export_dist2 = 1;
H.export_dist3 = 1;
H.export_dist4 = 1;
guidata(hObject,H);
plot_distribution(hObject, eventdata, H)

function autob_Callback(hObject, eventdata, H)
plot_distribution(hObject, eventdata, H)

function ydz_bins_Callback(hObject, eventdata, H)
set(H.autob,'Value',0)
plot_distribution(hObject, eventdata, H)

function sety_Callback(hObject, eventdata, H)
plot_distribution(hObject, eventdata, H)

function ymin_Callback(hObject, eventdata, H)
set(H.sety,'Value',0)
plot_distribution(hObject, eventdata, H)

function ymax_Callback(hObject, eventdata, H)
set(H.sety,'Value',0)
plot_distribution(hObject, eventdata, H)

function info_Callback(hObject, eventdata, H)

function copy_results_Callback(hObject, eventdata, H)
data = get(H.uitable6, 'Data');
copy(data);

function save_results_Callback(hObject, eventdata, H)

function cluster_sort_Callback(hObject, eventdata, H)
plot_distribution(hObject, eventdata, H)
