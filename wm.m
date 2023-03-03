function [WM, out1s, out2s, MSWD] = wm(dist_data)


t = sum(dist_data(:,1)./(dist_data(:,2).*dist_data(:,2))) / sum(1./(dist_data(:,2).*dist_data(:,2))); % Weighted Mean


		data2 = dist_data;
		data2(:,2) = data2(:,2).*2; % double the uncertainty to get the MSWD at 1 sigma.... THIS DOESN'T MAKE SENSE TO ME. SEEMS BACKWARDS!
		s = 1/sqrt(sum(1./(data2(:,2).*data2(:,2)))); % SE
		MSWD = 1/(length(data2(:,1))-1).*sum(((data2(:,1)- (sum(data2(:,1)./(data2(:,2).^2))/sum(1./(data2(:,2).^2))) ).^2)./((data2(:,2)./2).^2)); %MSWD at 1 sigma matches Isoplot


students_t = [12.71	4.303	3.182	2.776	2.571	2.447	2.365	2.306	2.262	2.228	2.201	2.179	2.16	2.145	2.131	2.12	2.11	2.101	2.093	2.086	2.08 ...
	2.074	2.069	2.064	2.06	2.056	2.052	2.048	2.045];

% 95% confidence interval using 2-sided Student's t
% if length(dist_data(:,1))-1 < 30
% 	conf95 = students_t(1,(length(dist_data(:,1))-1)) * s/2 *  sqrt(MSWD); 
% elseif length(dist_data(:,1))-1 >= 30 && length(dist_data(:,1))-1 < 40
% 	conf95 = 2.042 * s/2 *  sqrt(MSWD); 
% elseif length(dist_data(:,1))-1 >= 40 && length(dist_data(:,1))-1 < 50
% 	conf95 = 2.021 * s/2 *  sqrt(MSWD);
% elseif length(dist_data(:,1))-1 >= 50 && length(dist_data(:,1))-1 < 60
% 	conf95 = 2.009 * s/2 *  sqrt(MSWD);
% elseif length(dist_data(:,1))-1 >= 60 && length(dist_data(:,1))-1 < 80
% 	conf95 = 2.000 * s/2 *  sqrt(MSWD);
% elseif length(dist_data(:,1))-1 >= 80 && length(dist_data(:,1))-1 < 100
% 	conf95 = 1.99 * s/2 *  sqrt(MSWD);
% elseif length(dist_data(:,1))-1 >= 100 && length(dist_data(:,1))-1 < 120
% 	conf95 = 1.984 * s/2 *  sqrt(MSWD);
% elseif length(dist_data(:,1))-1 >= 120
% 	conf95 = 1.96 * s/2 *  sqrt(MSWD);
% end
% 
% y = conf95/sqrt(MSWD); %y at 2 sigma
% 
% z = y*sqrt(MSWD);

WM = t;
out1s= s/2;
out2s= s;
MSWD = MSWD;

