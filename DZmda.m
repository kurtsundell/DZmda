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
H.exportplot = 0;
guidata(hObject,H);

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

%%%% All MDA calculations assume 1 sigma % input below. Adjusted above if needed.

%keep original data
dist_data_perm = dist_data;
dist_data = sortrows(dist_data,1);

xint = 0.1;

x=0:xint:4500;

% Make PDP
pdp = pdp5(dist_data(:,1),dist_data(:,2),0,4500,xint); %1 sigma pdp input data
pdp = pdp/1/sum(pdp); % normalize pdp to 1

% MDA calculations
spcy_s = 0.0003; %adjust scatter triangles up

% Youngest single grain (YSG)
YSG_hi_tmp = dist_data(:,1) + dist_data(:,2);
[~,sortIdx0] = sort(YSG_hi_tmp);
dist_data_sort0 = dist_data(sortIdx0,:);
YSG = dist_data_sort0(1,1);
YSG_1s = dist_data_sort0(1,2);
YSG_2s = 2*dist_data_sort0(1,2);

% Youngest graphical peak of PDP (YPP)
[pks,locs] = findpeaks(pdp,x);
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

if get(H.setx,'Value') == 1
	xmin = round(YGF - 10*YGF_1s);
	xmax = round(YGF + 10*YGF_1s);
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
	[YSP_wm(i,1),YSP_1s(i,1),YSP_2s(i,1),YSP_mswdt(i,1)] = wm(YSP_data);
end
YSP_sub = 1 - YSP_mswdt;
[~,YSP_idx] = min(abs((YSP_sub)));
YSP = YSP_wm(YSP_idx,1);
YSP_1s = YSP_1s(YSP_idx,1);
YSP_2s = YSP_2s(YSP_idx,1);
YSP_mswd = YSP_mswdt(YSP_idx,1);
if min(YSP_mswdt) > 1
	clear YSP_sub
	ysptest = 0;
	while ysptest < 1
		for j = 2:length(dist_data(:,1))-1
			dist_data_tmp = dist_data(j:end,:);
			for i = 1:length(dist_data_tmp(:,1))-1
				dist_data_tmp_in = dist_data_tmp(1:i+1,:);
				[YSP_wmB(i,1),YSP_1sB(i,1),YSP_2sB(i,1),YSP_mswdtB(i,1)] = wm(dist_data_tmp_in);
			end
			YSP_sub = 1 - YSP_mswdtB;
			[~,YSP_idx] = min(abs((YSP_sub)));
			YSP = YSP_wmB(YSP_idx,1);
			YSP_1s = YSP_1sB(YSP_idx,1);
			YSP_2s = YSP_2sB(YSP_idx,1);
			YSP_mswd = YSP_mswdtB(YSP_idx,1);
			if min(YSP_mswdtB) < 1
				ysptest = 2;
				break
			end
			if j == length(dist_data(:,1))-1
				YSP = 0;
				YSP_1s = 0;
				YSP_2s = 0;
				YSP_mswd = 0;
				ysptest = 2;
				break
			end
			clear YSP_wmB YSP_1sB YSP_2sB YSP_mswdtB dist_data_tmp dist_data_tmp_in YSP_sub
		end
	end
end

% Make Age plot
if H.exportplot == 1
	figure;
end
if H.exportplot == 0
	cla(H.axes_comp,'reset');
	axes(H.axes_comp);
end
hold on

gauss_adj = Yhat*(1/(max(Yhat)/pks(1,1)));
p = plot(x, pdp, 'Color', 'b', 'LineWidth', 2);
p1 = plot(x32,gauss_adj,'Color',[1 0 0]);
set(p,'linewidth',2)
set(p1,'linewidth',2)
s1 = scatter(locs(1,1),pks(1,1), 150, 'filled', 'v', 'markeredgecolor', 'k',  'markerfacecolor', 'b', 'linewidth', 2);
s2 = scatter(YGF,pks(1,1), 150, 'filled', 'v', 'markeredgecolor', 'k',  'markerfacecolor', [1 0 0], 'linewidth', 2);
lgnd = legend([p, p1, s1, s2], 'Probability Density Plot', 'Youngest Gaussian Fit', 'YPP', 'YGF', 'fontsize',12);
set(lgnd,'Color','w');
set(gca, 'box', 'on')
set(gca, 'fontsize',14)
xlabel('Age (Ma)','Color','k')
ylabel('Probability','Color','k')
xlim([xmin xmax])

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

if H.exportplot == 1
	figure;
end
if H.exportplot == 0
	cla(H.axes58,'reset');
	axes(H.axes58);
end
hold on

ppp2 = plot([x2; x2], [(data2(:,1)+data2(:,2))'; (data2(:,1)-data2(:,2))'], '-r', 'Color', [.4 .6 1], 'LineWidth',5); % Error bars, nicer than the errorbar function
ppp1 = plot([x2; x2], [(data2(:,1)+data2(:,2)./2)'; (data2(:,1)-data2(:,2)./2)'], '-r', 'Color', 'k', 'LineWidth',5); % Error bars, nicer than the errorbar function
scatter(x2, data2(:,1), 150, 'filled','d', 'markeredgecolor','k','markerfacecolor','w','linewidth',1)
legend([ppp2(1) ppp1(1)], [{'2s'}, {'1s'}], 'Location','southwest');
xlim([xmin2 xmax2])
ylim([xmin xmax])

if H.exportplot == 0
	set(gca,'YTickLabel',[])
end

if H.exportplot == 1
	ylabel('Age (Ma)')
end

xlabel('Rank Number')
set(gca, 'box', 'on')
set(gca, 'fontsize',14)

camroll(-90) % rotate 90 deg
if H.exportplot == 1
	set(gca, 'YAxisLocation', 'right')
end

% Make results plot
if H.exportplot == 1
	figure;
end
if H.exportplot == 0
	cla(H.axes_mda,'reset');
	axes(H.axes_mda);
	
end
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
Y = [YSG;YPP;YGF;YGC1s;YGC2s;Y3Zo;Y3Za;74.1;75];
Y_unc = [YSG_2s;0.000001;YGF_2s;YGC1s_2s;YGC2s_2s;Y3Zo_2s;Y3Za_2s;0.3;1.1];

if get(H.sety,'Value') == 1
	yminp = round(min(nonzeros(Y(:,1))-nonzeros(Y_unc(:,1))) - 1); % make nice plots
	ymaxp = round(max(nonzeros(Y(:,1))+nonzeros(Y_unc(:,1))) + 1); % make nice plots
	set(H.ymin,'String',yminp)
	set(H.ymax,'String',ymaxp)
end

if get(H.sety,'Value') == 0
	yminp = str2num(get(H.ymin,'String'));
	ymaxp = str2num(get(H.ymax,'String'));
end

x3 = [1:1:9];
pp2 = plot([x3; x3], [(Y+Y_unc)'; (Y-Y_unc)'], '-r', 'Color', [.4 .6 1], 'LineWidth',7); % Error bars, nicer than the errorbar function
pp1 = plot([x3; x3], [(Y+Y_unc./2)'; (Y-Y_unc./2)'], '-r', 'Color', 'k', 'LineWidth',7); % Error bars, nicer than the errorbar function
scatter(X, Y, 250, 'filled','s', 'markeredgecolor','k','markerfacecolor','w','linewidth',2)
ylabel('Age (Ma)')
xlim([0 10])
ylim([yminp ymaxp])
set(gca, 'box', 'on')
names = {'YSG'; 'YPP'; 'YGF'; 'YGC1s'; 'YGC2s'; 'Y3Zo'; 'Y3Za'; 'TAU'; 'YSP'};
set(gca,'xtick',[1:9],'xticklabel',names)
set(gca, 'fontsize',14)
xlabel('Method')
set(gca,'YDir','reverse')
legend([pp2(1) pp1(1)], [{'2s'}, {'1s'}], 'Location','southwest');

% Make table
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
	tdata(4,1:3) = {'n<2'};
	tdata(4,5) = {'n<2'};
end
if length(YGC2s_data(:,1)) < 3
	tdata(5,1:3) = {'n<3'};
	tdata(5,5) = {'n<3'};
end
if Y3Zo == 0
	tdata(6,1:3) = {'n<3'};
	tdata(6,5) = {'n<3'};
end
if length(TAU_data(:,1)) < 3
	tdata(8,1:3) = {'n<3'};
	tdata(8,5) = {'n<3'};
end
if YSP == 0
	tdata(9,1:2) = {'NA'};
end

set(H.uitable6,'data',tdata)

H.exportplot = 0;
guidata(hObject,H);

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

function autob_Callback(hObject, eventdata, H)
plot_distribution(hObject, eventdata, H)

function sety_Callback(hObject, eventdata, H)
plot_distribution(hObject, eventdata, H)

function ymin_Callback(hObject, eventdata, H)
set(H.sety,'Value',0)
plot_distribution(hObject, eventdata, H)

function ymax_Callback(hObject, eventdata, H)
set(H.sety,'Value',0)
plot_distribution(hObject, eventdata, H)

function copy_results_Callback(hObject, eventdata, H)
data = get(H.uitable6, 'Data');
copy(data);

function save_results_Callback(hObject, eventdata, H)

function cluster_sort_Callback(hObject, eventdata, H)
plot_distribution(hObject, eventdata, H)

function exportplots_Callback(hObject, eventdata, H)
H.exportplot = 1;
guidata(hObject,H);
plot_distribution(hObject, eventdata, H)
