function [] = cylcomp_plotWF(name, waveG, waveC, waveF, wec, r0, As, varargin)

[opts, args] = checkOptions({{'clims', 3}, {'fontsize', 1}, {'figdims', 1}}, varargin);


scale = 100;

if (opts(1))
    clim1 = args{1}{1};
    clim2 = args{1}{2};
    clim3 = args{1}{3};
else
    clim1 = [0.7 1.3];
    clim2 = scale*[0 0.01];
    clim3 = scale*[0 0.1];
end

if (opts(2))
    fontsi = args{2};
else
    fontsi = 8;
end

if (opts(3))
    figdims = args{3};
    figwid = figdims(1);
    fighei = figdims(2);
else
    figwid = 6.5;
    fighei = 8.5;
end
    
tmar = 0.025;
lmar = 0.16;
lrspc = 0.07;
tbspc = 0.0;

pwid = 0.3;
phei = 0.3;

figure;
set(gcf, 'PaperPosition', [0 0 figwid fighei]);


[X, Y] = waveG.FieldPoints;

thetar = 0:2*pi/100:2*pi;
cirx = r0*cos(thetar);
ciry = r0*sin(thetar);

waveG.BodyMotions = wec.Motions;
waveC.BodyMotions = wec.Motions;
waveF.BodyMotions = wec.Motions;

% WAMIT
eta = waveG.Elevation('Total');
etaG = eta{1};

subplot('position', [lmar 1-tmar-phei pwid phei])
pcolor(X,Y,abs(etaG));
fet;
set(gca, 'clim', clim1, 'xticklabel', [], 'yticklabel', [], 'fontsize', fontsi)
ylabel('WAMIT', 'fontsize', fontsi)
title('Wave fields', 'fontsize', fontsi);

% WF colorbar
cbx = lmar + pwid + 0.04;
cby = 1-tmar-phei + 0.04;
chei = phei-0.08;
cwid = 0.04;

caxis = colorbar('location','manual','position',[cbx cby cwid chei]);
set(caxis, 'fontsize', fontsi);
%xlabel(caxis, '|\eta/a|', 'fontsize', fontsi);
x1 = caxis.Position(1);
y1 = caxis.Position(2);
dxcl = -0.06;
dycl = -0.08;
annotation(gcf, 'textbox', [x1+dxcl+0.02 y1+dycl+0.01 0.2 0.07], 'string', '|\eta/a|', 'linestyle', 'none', 'fontsize', fontsi);

% Cylinder
eta = waveC.Elevation('Total');
etaC = eta{1};

subplot('position', [lmar 1-tmar-2*phei-tbspc pwid phei])
pcolor(X,Y,abs(etaC));
fet;
set(gca, 'clim', clim1, 'xticklabel', [], 'yticklabel', [], 'fontsize', fontsi)
ylabel('Cylindrical', 'fontsize', fontsi)

hold on;
plot(cirx, ciry, 'w');

subplot('position', [lmar+pwid+lrspc 1-tmar-2*phei-tbspc pwid phei])
pcolor(X,Y,scale*abs(etaC-etaG)./abs(etaG));
%pcolor(X,Y,abs(etaC-etaG));
fet;
set(gca, 'clim', clim2, 'xticklabel', [], 'yticklabel', [], 'fontsize', fontsi)
title('Error', 'fontsize', fontsi);

hold on;
plot(cirx, ciry, 'w');

% Error colorbar 1
cbx = lmar + 2*pwid + lrspc + 0.04;
cby = 1-tmar-2*phei-tbspc + 0.04;

caxis = colorbar('location','manual','position',[cbx cby cwid chei]);
set(caxis, 'fontsize', fontsi);
%xlabel(caxis, 'Err (%)', 'fontsize', fontsi);
x1 = caxis.Position(1);
y1 = caxis.Position(2);
annotation(gcf, 'textbox', [x1+dxcl y1+dycl 0.2 0.07], 'string', 'Err (%)', 'linestyle', 'none', 'fontsize', fontsi);

% Far-field
eta = waveF.Elevation('Total');
etaF = eta{1};

subplot('position', [lmar 1-tmar-3*phei-2*tbspc pwid phei])
pcolor(X,Y,abs(etaF));
fet;
set(gca, 'clim', clim1, 'xtick', [-20 0 20], 'ytick', [-20 0 20], 'fontsize', fontsi)
ylabel(['Far-field (y/d)'], 'fontsize', fontsi)
xlabel('(x/d)', 'fontsize', fontsi);

hold on;
plot(cirx, ciry, 'w');

subplot('position', [lmar+pwid+lrspc 1-tmar-3*phei-2*tbspc pwid phei])
pcolor(X,Y,scale*abs(etaF-etaG)./abs(etaG));
%pcolor(X,Y,abs(etaF-etaG)./abs(etaG));
fet;
set(gca, 'clim', clim3, 'xticklabel', [], 'yticklabel', [], 'fontsize', fontsi)

hold on;
plot(cirx, ciry, 'w');

% Error colorbar 2
cbx = lmar + 2*pwid + lrspc + 0.04;
cby = 1-tmar-3*phei-2*tbspc + 0.04;

caxis = colorbar('location','manual','position',[cbx cby cwid chei]);
set(caxis, 'fontsize', fontsi);
%xlabel(caxis, 'Err (%)', 'fontsize', fontsi);
x1 = caxis.Position(1);
y1 = caxis.Position(2);
annotation(gcf, 'textbox', [x1+dxcl y1+dycl 0.2 0.07], 'string', 'Err (%)', 'linestyle', 'none', 'fontsize', fontsi);

% colormap
cmap = colormap(gca);
ncol = size(cmap,1);
cstep = floor(ncol/10);
imap = 1:cstep:ncol;
cmap2 = cmap(imap,:);
colormap(cmap2)

% M and N text
[M N] = size(As{1});
M = M - 1;
Ns = (N - 1)/2;
N = size(As{2},2);
Nr = (N - 1)/2;

if (length(As) > 2)
    N = size(As{3},2);
    Nr2 = (N - 1)/2;
else
    Nr2 = -1;
end

if (Nr2 < 0)
    annotation(gcf, 'textbox', [0.7 0.89 0.1 0.1], 'string', [{['L = ' num2str(M)]}, {' '}, {['M^S = ' num2str(Ns)]}, {' '}, {['M^R = ' num2str(Nr)]}], 'linestyle', 'none', 'fontsize', fontsi)
else
   annotation(gcf, 'textbox', [0.7 0.89 0.1 0.1], 'string', [{['L = ' num2str(M)]}, {' '}, {['M^S = ' num2str(Ns)]}, {' '}, {['M^{R}_{1} = ' num2str(Nr)]}, {' '}, {['M^{R}_{2} = ' num2str(Nr2)]}], 'linestyle', 'none', 'fontsize', fontsi)
end


print('-dpng', '-r400', name);


end

