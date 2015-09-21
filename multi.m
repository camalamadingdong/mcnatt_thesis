%% Multi-body interactions
% copied from jfm_PaperPics in N:\RDS\WAMIT\array_study\jfm
% 0) set up stuff

folder = 'C:\Users\s1213969\Dropbox\matlab\thesis\';
fontsi = 10;
figwid = 14;

%% 1) pics of geometries

[~, cyl] = multi_makeRun('Cyl', 'Geo', [], [] , []);

[~, atten] = multi_makeRun('Atten', 'Geo', [], [] , []);

fighei = 5;
fontsi = 9;

figure; 
set(gcf, 'PaperPosition', [0 0 figwid fighei]);

phei1 = 0.8;
pwid1 = 0.15;

phei2 = phei1;
pwid2 = 0.5;

lrspc = 0.2;
lmar = 0.09;
tmar = 0.01;

figure;
set(gcf, 'PaperPosition', [0 0 figwid fighei]);


subplot('Position', [lmar 1-tmar-phei1 pwid1 phei1]);

plot(cyl.PanelGeo, 'OnlyWet');
set(gca, 'view', [-37.5000 30]);
axis equal
set(gca, 'FontSize', fontsi, 'xtick', [-1 0 1], 'ytick', [-1 0 1], 'ztick', [-0.5 0]);
%ylabel('y/a', 'FontSize', fontsi);
annotation(gcf,'textbox',[0.23 0.36 0.2 0.06],'String',{'x/d'}, 'fontsize', fontsi,'LineStyle','none');
annotation(gcf,'textbox',[0.03 0.4 0.2 0.06],'String',{'y/d'}, 'fontsize', fontsi,'LineStyle','none');
zlabel('z/d', 'FontSize', fontsi);
set(get(gca, 'zlabel' ), 'Rotation' ,0 )

annotation(gcf,'textbox',[0.1 0.12 0.2 0.06],'String',{'a) Cylinder'}, 'fontsize', fontsi,'LineStyle','none');

subplot('Position', [lmar+pwid1+lrspc 1-tmar-phei2 pwid2 phei2]);

plot(atten.PanelGeo, 'OnlyWet');
plot3([atten.HingePos(1) atten.HingePos(1)], [-0.5 0.5], [0 0], 'k', 'linewidth', 2);
plot3([atten.HingePos(2) atten.HingePos(2)], [-0.5 0.5], [0 0], 'k', 'linewidth', 2);
axis equal;
set(gca, 'FontSize', fontsi, 'xtick', [-5 0 5], 'ytick', [-0.5 0.5], 'ztick', [-0.4 0]);
annotation(gcf,'textbox',[0.7 0.36 0.2 0.06],'String',{'x/d'}, 'fontsize', fontsi,'LineStyle','none');
annotation(gcf,'textbox',[0.32 0.4 0.2 0.06],'String',{'y/d'}, 'fontsize', fontsi,'LineStyle','none');
zlabel('z/d', 'FontSize', fontsi);
set(get(gca, 'zlabel' ), 'Rotation' ,0 )
view(gca,[-21.5 18]);
view(gca,[-25 18]);

annotation(gcf,'textbox',[0.62 0.12 0.2 0.06],'String',{'b) Attenuator'}, 'fontsize', fontsi,'LineStyle','none');

print('-dpng', '-r400', [folder 'pics\geos']);

%% 2) performance pic

lam = 0.4:0.2:20;

beta = 0;

cylComp = multi_makeRun('Cyl', 'Perf', [folder 'runs'], lam, beta);

attenComp = multi_makeRun('Atten', 'Perf', [folder 'runs'], lam, beta);

fontsi = 10;

Xi = squeeze(cylComp.Motions);

Xihc = squeeze(Xi(:,3));
Xisc = squeeze(Xi(:,1));
Xipc = squeeze(Xi(:,5));

D = zeros(8,8);
attenComp.SetDpto(D);
Xi = squeeze(attenComp.Motions);

Xiha = squeeze(Xi(:,3));
Xisa = squeeze(Xi(:,1));
Xipa = squeeze(Xi(:,5));
Xif1 = squeeze(Xi(:,7));
Xif2 = squeeze(Xi(:,8));

Ef = IWaves.UnitEnergyFlux(1000, attenComp.T, attenComp.H).';

D(7,7) = 35000;
D(8,8) = 35000;
attenComp.SetDpto(D);
pow = squeeze(attenComp.Power);
CWa = squeeze(pow(:,7)+pow(:,8))./Ef.';

nondimD = 35000/(1000*sqrt(1)*sqrt(IWaves.G));

Xi = squeeze(attenComp.Motions);

Xihad = squeeze(Xi(:,3));
Xisad = squeeze(Xi(:,1));
Xipad = squeeze(Xi(:,5));
Xif1d = squeeze(Xi(:,7));
Xif2d = squeeze(Xi(:,8));


fighei = 10;

figure;
set(gcf, 'PaperPosition', [0 0 figwid fighei]);

phei1 = 0.30;
pwid1 = 0.29;

phei2 = phei1;
pwid2 = pwid1;

lmar = 0.09;
tmar = 0.15;
lrspc = 0.23;
tbspc1 = 0.13;

subplot('Position', [lmar 1-tmar-phei1 pwid1 phei1]);

[Ax, h1, h2] = plotyy(lam,abs([Xisc, Xihc]),lam,180/pi*abs(Xipc));
set(Ax(1), 'ytick', [0 2 4 6], 'FontSize', fontsi);
set(Ax(2), 'ytick', [0 1 2 3], 'ycolor', 'k', 'FontSize', fontsi);
set(get(Ax(1),'ylabel'),'String',{'Surge, Heave (m/m)'},'FontSize', fontsi); 
set(get(Ax(2),'ylabel'),'Color', 'k','String',{'Pitch (deg/m)'},'Color', 'k','FontSize', fontsi,'Color',[1 0 0]);
set(get(Ax(2),'ylabel'), 'Color', [0 0 0]);
xlabel('\lambda/d','FontSize', fontsi);
title('Cylinder', 'FontSize', fontsi);

annotation(gcf,'textbox',[0.4 0.93 0.3 0.06],'String',{'Undamped'}, 'fontsize', fontsi,'LineStyle','none','FontAngle','italic');

subplot('Position', [lmar+pwid1+lrspc 1-tmar-phei1 pwid1 phei1]);

[Ax, h1, h2] = plotyy(lam,abs([Xisa, Xiha]),lam,180/pi*abs([Xipa, Xif1, Xif2]));
set(Ax(1), 'xticklabel', '', 'FontSize', fontsi);
set(Ax(2), 'xticklabel', '','FontSize', fontsi);
set(get(Ax(1),'Ylabel'),'String',{'Surge, Heave (m/m)'},'FontSize', fontsi) 
set(get(Ax(2),'Ylabel'),'String',{'Pitch, Flex (deg/m)'},'FontSize', fontsi) 
title('Attenuator', 'FontSize', fontsi);

legax = legend({'Surge', 'Heave', 'Pitch', 'Forward Flex', 'Aft Flex'},'FontSize', fontsi);
set(legax,'Position', [.16,.16,.1,.2]);

subplot('Position', [lmar+pwid1+lrspc 1-tmar-2*phei1-tbspc1 pwid1 phei1]);

[Ax, h1, h2] = plotyy(lam,abs([Xisad, Xihad]),lam,180/pi*abs([Xipad, Xif1d, Xif2d]));
set(Ax(1), 'FontSize', fontsi);
set(Ax(2), 'FontSize', fontsi);
set(get(Ax(1),'ylabel'),'String',{'Surge, Heave (m/m)'}, 'fontsize', fontsi) 
set(get(Ax(2),'Ylabel'),'String',{'Pitch, Flex (deg/m)'},'fontsize', fontsi) 
xlabel('\lambda/d','FontSize', fontsi);

annotation(gcf,'textbox',[0.42 0.45 0.3 0.06],'String',{'Damped'}, 'fontsize', fontsi,'LineStyle','none','FontAngle','italic');


print('-dpng', '-r400', [folder 'pics\motions_fig']);

%% 3) compute hydro bodies

lam = [3, 10, 30];

[~, ~, cylHydBod] = multi_makeRun('Cyl', 'HydroBody', [folder 'runs'], lam, []);

save([folder '\multi\cyl_hydBod'], 'cylHydBod');

[~, ~, attenHydBod] = multi_makeRun('Atten', 'HydroBody', [folder 'runs'], lam, []);

save([folder '\multi\atten_hydBod'], 'attenHydBod');

%% 4.1) WAMIT distance calcs

lam = [3, 10, 30];
beta = [0 pi/4 pi/2];

% Cylinder
dis1 = 2 + 0.4;
delta = 0.2;
dis2 = 20;
Dist = dis1:delta:dis2;

Aw = zeros(length(Dist), length(lam), 12, 12);
Bw = zeros(length(Dist), length(lam), 12, 12);
Fexw = zeros(length(Dist), length(lam), length(beta), 12);
Xiw = zeros(length(Dist), length(lam), length(beta), 12);

for n = 1:length(Dist)
    disp(['******************************************** D = ' num2str(Dist(n)) ' ********************************************'])

    comp = multi_makeRun('Cyl', 'Dist', [folder 'runs'], lam, beta, Dist(n), 'Side2Side');
    
    Aw(n,:,:,:) = squeeze(comp.A);
    Bw(n,:,:,:) = squeeze(comp.B);
    Fexw(n,:,:,:) = squeeze(comp.Fex);
    Xiw(n,:,:,:) = comp.Motions;
end

save([folder '\multi\cyl_dist'], 'Dist', 'Aw', 'Bw', 'Fexw', 'Xiw');

% Attenuator, side to side
dis1 = 6.2;
delta = 0.2;
dis2 = 20;
Dist = dis1:delta:dis2;

Aw = zeros(length(Dist), length(lam), 16, 16);
Bw = zeros(length(Dist), length(lam), 16, 16);
Fexw = zeros(length(Dist), length(lam), length(beta), 16);
Xiw = zeros(length(Dist), length(lam), length(beta), 16);

for n = 1:length(Dist)
    disp(['******************************************** D = ' num2str(Dist(n)) ' ********************************************'])

    comp = multi_makeRun('Atten', 'Dist', [folder 'runs'], lam, beta, Dist(n), 'Side2Side');
    
    Aw(n,:,:,:) = squeeze(comp.A);
    Bw(n,:,:,:) = squeeze(comp.B);
    Fexw(n,:,:,:) = squeeze(comp.Fex);
    Xiw(n,:,:,:) = comp.Motions;
end

save([folder '\multi\atten_distSide'], 'Dist', 'Aw', 'Bw', 'Fexw', 'Xiw');

% Attenuator, front to back
dis1 = 10.2;
delta = 0.2;
dis2 = 20;
Dist = dis1:delta:dis2;

Aw = zeros(length(Dist), length(lam), 16, 16);
Bw = zeros(length(Dist), length(lam), 16, 16);
Fexw = zeros(length(Dist), length(lam), length(beta), 16);
Xiw = zeros(length(Dist), length(lam), length(beta), 16);

for n = 1:length(Dist)
    disp(['******************************************** D = ' num2str(Dist(n)) ' ********************************************'])

    comp = multi_makeRun('Atten', 'Dist', [folder 'runs'], lam, beta, Dist(n), 'Front2Back');
    
    Aw(n,:,:,:) = squeeze(comp.A);
    Bw(n,:,:,:) = squeeze(comp.B);
    Fexw(n,:,:,:) = squeeze(comp.Fex);
    Xiw(n,:,:,:) = comp.Motions;
end

save([folder '\multi\atten_distFront'], 'Dist', 'Aw', 'Bw', 'Fexw', 'Xiw');

%% 4.2) IT Cylinder distance

load cyl_dist

load cyl_hydBod

hydBod = cylHydBod;
side2side = true;
skipDist = 0;

beta = [0 pi/4 pi/2];

[Aa, Ba, Fexa, Xia] = multi_computeHBdist(hydBod, Dist, beta, side2side, skipDist);



% Added Mass - A28, A39, A410
% Fex - Beta = 0: Fex2, Fex3, Fex4

%% 4.3) Cylinder Distance plot


fighei = 16;

rho = 1000;
a = 1;
g = IWaves.G;

dimA = 1./[rho*a^3, rho*a];
dimF = 1./[rho*a^3*g, rho*a^4*g];

figure
set(gcf, 'PaperPosition', [0 0 figwid fighei]);

col1 = [0 0.4470 0.7410];
col2 = [0.85 0.325 0.098];
col3 = [0.929 0.694 0.125];

%--------------------------------------------------------------------------
Nrow = 3;
M1 = [2 8; 3 9; 4 10];
DimI = [1, 1, 2];

pwid = 0.35;
phei = 1/Nrow - 0.125;
lmar = 0.13;
lrspc = 0.13;
tbspc1 = 0.04;
tmar = 0.05;

lwida = 0.5;
lwidw = 0.8;

for rowI = 1:Nrow
    m = M1(rowI,1);
    n = M1(rowI,2);
    dimI = DimI(rowI);
    
    subplot('Position', [lmar 1-tmar-(rowI-1)*(phei+tbspc1)-phei pwid phei]);
    %subplot(Nrow,2,(rowI-1)*2+1);
    
    plot(Dist, dimA(dimI)*squeeze(Aw(:,1,m,n)), 'LineStyle','--', 'Color', col1);
    hold on;
    plot(Dist, dimA(dimI)*squeeze(Aa(:,1,m,n)), 'Color', col1);
    plot(Dist, dimA(dimI)*squeeze(Aw(:,2,m,n)),'LineStyle','--', 'Color', col2);
    plot(Dist, dimA(dimI)*squeeze(Aa(:,2,m,n)), 'Color', col2);
    plot(Dist, dimA(dimI)*squeeze(Aw(:,3,m,n)), 'LineStyle','--', 'Color', col3);
    plot(Dist, dimA(dimI)*squeeze(Aa(:,3,m,n)), 'Color', col3);
    set(gca, 'xlim', [2 10],'xtick',[2 4 6 8 10], 'fontsize', fontsi);
    
    if (rowI == 1)
        title(['Added Mass'], 'fontsize', fontsi);
    end
    
    if (rowI == Nrow)    
        xlabel('Distance (\delta/d)','fontsize', fontsi);
        legax = legend({'WAM, \lambda/d=3', 'IT', 'WAM, \lambda/d=10', 'IT', 'WAM, \lambda/d=30', 'IT'},'FontSize', fontsi);
        set(legax,'Position', [0.5, 0.0, 0.1, 0.2],'Box', 'off');
    else
        set(gca, 'xticklabel', '');
    end
    
    if (rowI == 3)
        ylabel(['A_{' num2str(m) ',' num2str(n) '} /\rho d'],'fontsize', fontsi);
    else
        ylabel(['A_{' num2str(m) ',' num2str(n) '} /\rho d^{3}'],'fontsize', fontsi);
    end
    
    subplot('Position', [lmar+pwid+lrspc 1-tmar-(rowI-1)*(phei+tbspc1)-phei pwid phei]);
    m = 1;
    plot(Dist, dimF(dimI)*abs(squeeze(Fexw(:,1,m,n))), 'LineStyle','--', 'Color', col1);
    hold on;
    plot(Dist, dimF(dimI)*abs(squeeze(Fexa(:,1,m,n))), 'Color', col1);
    plot(Dist, dimF(dimI)*abs(squeeze(Fexw(:,2,m,n))), 'LineStyle','--', 'Color', col2);
    plot(Dist, dimF(dimI)*abs(squeeze(Fexa(:,2,m,n))), 'Color', col2);
    plot(Dist, dimF(dimI)*abs(squeeze(Fexw(:,3,m,n))), 'LineStyle','--', 'Color', col3);
    plot(Dist, dimF(dimI)*abs(squeeze(Fexa(:,3,m,n))), 'Color', col3);
    set(gca, 'xlim', [2 10], 'xtick',[2 4 6 8 10], 'fontsize', fontsi);
    
    if (rowI == 1)
        title(['Force (\beta = 0)'], 'fontsize', fontsi);
    end
    if (rowI == 3)
        ylabel(['f_{' num2str(rowI+1) '} / \rho g d^4'],'fontsize', fontsi);
    else
        ylabel(['f_{' num2str(rowI+1) '} / \rho g d^3'],'fontsize', fontsi);
    end
    if (rowI == Nrow)    
        xlabel('Distance (\delta/d)','fontsize', fontsi);
    else
        set(gca, 'xticklabel', '');
    end
end

print('-dpng', '-r400', [folder 'pics\cyl_dist']);

%% 4.4) IT Atten distance

load atten_distSide

load atten_hydBod

beta = [0 pi/4 pi/2];

skipDist = 5.5;

Dist1 = Dist;
Aw1 = Aw;
Bw1 = Bw;
Fexw1 = Fexw;
Xiw1 = Xiw;

[Aa1, Ba1, Fexa1, Xia1] = multi_computeHBdist(attenHydBod, Dist, beta, true, skipDist);

% Added Mass - A22, A315
% Fex - Beta = pi/4: Fex2, Fex7
load atten_distFront

load atten_hydBod

beta = [0 pi/4 pi/2];

skipDist = 10;

Dist2 = Dist;
Aw2 = Aw;
Bw2 = Bw;
Fexw2 = Fexw;
Xiw2 = Xiw;

[Aa2, Ba2, Fexa2, Xia2] = multi_computeHBdist(attenHydBod, Dist, beta, false, skipDist);

% Added Mass - A19
% Fex - Beta = pi/2: Fex5

%% 4.5) Atten Distance plot

fighei = 18;

rho = 1000;
a = 1;
g = IWaves.G;

dimA = 1./[rho*a^3, rho*a^2, rho*a];
dimF = 1./[rho*a^3*g, rho*a^4*g];

figure
set(gcf, 'PaperPosition', [0 0 figwid fighei]);

%-------------------------------------------------------------------------

col1 = [0 0.4470 0.7410];
col2 = [0.85 0.325 0.098];
col3 = [0.929 0.694 0.125];

%pwid = 0.24;
pwid = 0.34;
phei = 0.18;
lmar = 0.14;
lrspc = 0.13;
tbspc1 = 0.1;
tbspc2 = 0.04;
tmar = 0.08;

% Plot 3
subplot('Position', [lmar 1-tmar-phei pwid phei]);

m = 1;
n = 9;
dimI = 1;

plot(Dist2, dimA(dimI)*squeeze(Aw2(:,1,m,n)), 'LineStyle','--', 'Color', col1);
hold on;
plot(Dist2, dimA(dimI)*squeeze(Aa2(:,1,m,n)), 'Color', col1);
plot(Dist2, dimA(dimI)*squeeze(Aw2(:,2,m,n)), 'LineStyle','--', 'Color', col2);
plot(Dist2, dimA(dimI)*squeeze(Aa2(:,2,m,n)), 'Color', col2);
plot(Dist2, dimA(dimI)*squeeze(Aw2(:,3,m,n)), 'LineStyle','--', 'Color', col3);
plot(Dist2, dimA(dimI)*squeeze(Aa2(:,3,m,n)), 'Color', col3);

set(gca, 'xlim', [10 20],'xtick',[10 15 20], 'fontsize', fontsi);
title(['Added Mass'],'fontsize', fontsi);
ylabel(['A_{' num2str(m) ',' num2str(n) '} /\rho d^{3}'],'fontsize', fontsi);

% Plot 2
subplot('Position', [lmar+pwid+lrspc 1-tmar-phei pwid phei]);

m = 3;
n = 5;
dimI = 2;

plot(Dist2, dimF(dimI)*abs(squeeze(Fexw2(:,1,m,n))), 'LineStyle', '--', 'Color', col1);
hold on;
plot(Dist2, dimF(dimI)*abs(squeeze(Fexa2(:,1,m,n))), 'Color', col1);
plot(Dist2, dimF(dimI)*abs(squeeze(Fexw2(:,2,m,n))), 'LineStyle','--', 'Color', col2);
plot(Dist2, dimF(dimI)*abs(squeeze(Fexa2(:,2,m,n))), 'Color', col2);
plot(Dist2, dimF(dimI)*abs(squeeze(Fexw2(:,3,m,n))), 'LineStyle','--', 'Color', col3);
plot(Dist2, dimF(dimI)*abs(squeeze(Fexa2(:,3,m,n))), 'Color', col3);

set(gca, 'xlim', [10 20], 'xtick',[10 15 20], 'fontsize', fontsi);
title(['Force (\beta = \pi/2)'], 'fontsize', fontsi);
ylabel(['f_{5} / \rho g d^4'], 'fontsize', fontsi);

% Plot 3
subplot('Position',  [lmar 1-tmar-(phei+tbspc1)-phei pwid phei]);

m = 2;
n = 2;
dimI = 1;

plot(Dist1, dimA(dimI)*squeeze(Aw1(:,1,m,n)),'LineStyle', '--', 'Color', col1);
hold on;
plot(Dist1, dimA(dimI)*squeeze(Aa1(:,1,m,n)), 'Color', col1);
plot(Dist1, dimA(dimI)*squeeze(Aw1(:,2,m,n)), 'LineStyle','--', 'Color', col2);
plot(Dist1, dimA(dimI)*squeeze(Aa1(:,2,m,n)), 'Color', col2);
plot(Dist1, dimA(dimI)*squeeze(Aw1(:,3,m,n)), 'LineStyle','--', 'Color', col3);
plot(Dist1, dimA(dimI)*squeeze(Aa1(:,3,m,n)), 'Color', col3);

set(gca, 'xlim', [5 20],'xtick',[5 10 15 20],'xticklabel', {},'fontsize', fontsi);
title(['Added Mass'], 'fontsize', fontsi);
ylabel(['A_{' num2str(m) ',' num2str(n) '} /\rho d^{3}'],'fontsize', fontsi);

annotation(gcf,'textbox',[0.44 0.93 0.2 0.06],'String',{'Front to Back'}, 'fontsize', fontsi,'LineStyle','none','FontAngle','italic');

% Plot 4
subplot('Position', [lmar+pwid+lrspc 1-tmar-(phei+tbspc1)-phei pwid phei]);

m = 2;
n = 2;
dimI = 1;

plot(Dist1, dimF(dimI)*abs(squeeze(Fexw1(:,1,m,n))),'LineStyle','--', 'Color', col1);
hold on;
plot(Dist1, dimF(dimI)*abs(squeeze(Fexa1(:,1,m,n))), 'Color', col1);
plot(Dist1, dimF(dimI)*abs(squeeze(Fexw1(:,2,m,n))), 'LineStyle','--', 'Color', col2);
plot(Dist1, dimF(dimI)*abs(squeeze(Fexa1(:,2,m,n))), 'Color', col2);
plot(Dist1, dimF(dimI)*abs(squeeze(Fexw1(:,3,m,n))), 'LineStyle','--', 'Color', col3);
plot(Dist1, dimF(dimI)*abs(squeeze(Fexa1(:,3,m,n))), 'Color', col3);

set(gca, 'xlim', [5 20], 'xtick',[5 10 15 20], 'xticklabel', {}, 'fontsize', fontsi);
title(['Force (\beta = \pi/4)'], 'fontsize', fontsi);
ylabel(['f_{' num2str(n) '} / \rho g d^3'], 'fontsize', fontsi);

annotation(gcf,'textbox',[0.44 0.64 0.2 0.06],'String',{'Side by Side'}, 'fontsize', fontsi,'LineStyle','none','FontAngle','italic');

% Plot 5
subplot('Position',  [lmar 1-tmar-3*phei-tbspc1-tbspc2 pwid phei]);

m = 3;
n = 15;
dimI = 2;

plot(Dist1, dimA(dimI)*squeeze(Aw1(:,1,m,n)),'LineStyle','--', 'Color', col1);
hold on;
plot(Dist1, dimA(dimI)*squeeze(Aa1(:,1,m,n)), 'Color', col1);
plot(Dist1, dimA(dimI)*squeeze(Aw1(:,2,m,n)), 'LineStyle','--', 'Color', col2);
plot(Dist1, dimA(dimI)*squeeze(Aa1(:,2,m,n)), 'Color', col2);
plot(Dist1, dimA(dimI)*squeeze(Aw1(:,3,m,n)), 'LineStyle','--', 'Color', col3);
plot(Dist1, dimA(dimI)*squeeze(Aa1(:,3,m,n)), 'Color', col3);

set(gca, 'xlim', [5 20],'xtick',[5 10 15 20],'fontsize', fontsi);
xlabel('Distance (\delta/d)','fontsize', fontsi);
ylabel(['A_{' num2str(m) ',' num2str(n) '} /\rho d^{2}'],'fontsize', fontsi);

% Plot 6
subplot('Position', [lmar+pwid+lrspc 1-tmar-3*phei-tbspc1-tbspc2 pwid phei]);

m = 2;
n = 7;
dimI = 2;

plot(Dist1, dimF(dimI)*abs(squeeze(Fexw1(:,1,m,n))),'LineStyle','--', 'Color', col1);
hold on;
plot(Dist1, dimF(dimI)*abs(squeeze(Fexa1(:,1,m,n))), 'Color', col1);
plot(Dist1, dimF(dimI)*abs(squeeze(Fexw1(:,2,m,n))), 'LineStyle','--', 'Color', col2);
plot(Dist1, dimF(dimI)*abs(squeeze(Fexa1(:,2,m,n))), 'Color', col2);
plot(Dist1, dimF(dimI)*abs(squeeze(Fexw1(:,3,m,n))), 'LineStyle','--', 'Color', col3);
plot(Dist1, dimF(dimI)*abs(squeeze(Fexa1(:,3,m,n))), 'Color', col3);

set(gca, 'xlim', [5 20], 'xtick',[5 10 15 20], 'fontsize', fontsi);
xlabel('Distance (\delta/d)','fontsize', fontsi);
ylabel(['f_{' num2str(n) '} / \rho g d^4'],'fontsize', fontsi);

legax = legend({'WAM, \lambda/d=3', 'IT', 'WAM, \lambda/d=10', 'IT', 'WAM, \lambda/d=30', 'IT'},'FontSize', fontsi);
set(legax,'Position', [0.5, 0.0, 0.1, 0.2],'Box', 'off');

print('-dpng', '-r400', [folder 'pics\atten_dist']);

%% 5.1) Medium Array - Cylinder - WAMIT

lam = [3, 10, 30];
beta = [0, pi/4];

pos = zeros(16,2);
ind = 1;
for m = 1:4
    for n = 1:4
        pos(ind, :) = [-7.5+(m-1)*5, -7.5+(n-1)*5];
        ind = ind + 1;
    end
end

fieldArray = BemFieldArray([-30 -30 0], [0.4 0.4 1], [301 151 1]);

[wcomp, ~, ~, wWF] = multi_makeRun('Cyl', 'MedArray', [folder 'runs'], lam, beta, pos, fieldArray);

save([folder '\multi\cyl_medArray'], 'wcomp', 'wWF');

%% 5.2) Medium Array - Cylinder - IT
load 'cyl_medArray';


floatBs = wcomp.Bodies;
wWF.BodyMotions = wcomp.Motions;

load 'cyl_hydBod';

for n = 1:length(floatBs)
    hbs(n) = HydroBody(cylHydBod);
    hbs(n).XYpos = floatBs(n).XYpos;
end

acomps = HydroArrayComp(hbs, wcomp.IncWaves);

[X, Y] = wWF.FieldPoints;
tic 
aWFs = acomps.WaveField(true, X, Y, 'NoVel');
toc

aWFs.BodyMotions = acomps.Motions;

%% 5.3) Medium Array - Cylinder - Plot

%load cyl_med_array

wWF.BodyMotions = wcomp.Motions;
aWFs.BodyMotions = acomps.Motions;

hbs = acomps.Bodies;
[X, Y] = wWF.FieldPoints;

fighei = 10;

type = 'Total';

clims = [0.2 1.8; 0.2 1.8; 0.2 1.8];
clims2 = [0 2; 0 2; 0 2];

etaW = wWF.Elevation(type);
etaA = aWFs.Elevation(type);

figure
set(gcf, 'PaperPosition', [0 0 figwid fighei]);

pwid = 0.22;
phei = 0.21;
lmar = 0.125;
lrspc = 0.015;
tbspc1 = 0.02;
tmar = 0.14;

chei = phei - 0.04;
cwid = 0.02;
clrspc = 0.02;
cbotspc = 0.02;

nb = 1;

axiz = zeros(9,1);
iax = 1;
caxiz = zeros(3,2);

for m = 1:3 
    
    nt = m;
    bpos = 1-tmar-(m-1)*(phei+tbspc1)-phei;
        
    if (m == 3)
        nb = 2;
        nt = 2;
        bpos = 1-tmar-(m-1)*(phei+2.5*tbspc1)-phei;
    end
        
    subplot('Position', [lmar bpos pwid phei]);
    axiz(iax) = gca; 
    iax = iax + 1;
    pcolor(X, Y, abs(etaW{nt, nb}));
    multi_drawBodiesOnPlot(hbs)
    fet;
    shading interp
    set(gca, 'clim', clims(m,:), 'fontsize', fontsi);
    if (m == 1)
        title ('WAMIT', 'fontsize', fontsi);
        ylabel('\lambda/d = 3', 'fontsize', fontsi);
        set(gca, 'xticklabel', '', 'yticklabel', '');
    elseif (m == 2)
        ylabel('\lambda/d = 10', 'fontsize', fontsi);
        set(gca, 'xticklabel', '', 'yticklabel', '');
    else
        xlabel('x/d', 'fontsize', fontsi);
        ylabel({'\lambda/d = 10','y/d'}, 'fontsize', fontsi);
        set(gca, 'xtick', [0 40 80]);
    end
            
    subplot('Position', [lmar+pwid+lrspc bpos pwid phei]);
    axiz(iax) = gca; 
    iax = iax + 1;
    pcolor(X, Y, abs(etaA{nt, nb}));
    multi_drawBodiesOnPlot(hbs)
    fet;
    shading interp
    set(gca, 'clim', clims(m,:), 'fontsize', fontsi);
    set(gca, 'xticklabel', '', 'yticklabel', '');
    if (m == 1)
        title ('IT', 'fontsize', fontsi);
    end
    
    cbx = lmar+2*pwid+lrspc+clrspc;
    cby = bpos + cbotspc;

    caxis = colorbar('location','manual','position',[cbx cby cwid chei]);
    caxiz(m,1) = caxis;
    set(caxis, 'fontsize', fontsi);
    
    if (m == 3)        
        set(caxis, 'fontsize', fontsi);
        dxcl = -0.02;
        dycl = -0.08;
        x1 = caxis.Position(1);
        y1 = caxis.Position(2);
        annotation(gcf, 'textbox', [x1+dxcl y1+dycl 0.2 0.07], 'string', '|\eta/a|', 'linestyle', 'none', 'fontsize', fontsi);
        %xlabel(caxis, '|\eta/a|', 'fontsize', fontsi);
    end
    
    subplot('Position', [lmar+2*pwid+7*lrspc bpos pwid phei]);
    axiz(iax) = gca; 
    iax = iax + 1;
    %pcolor(X, Y, abs((etaA{m,nb}-etaW{m,nb})./etaW{m,nb}));
    pcolor(X, Y, 100*abs((etaA{nt,nb}-etaW{nt,nb})));
    multi_drawBodiesOnPlot(hbs)
    fet;
    shading interp
    set(gca, 'clim', clims2(m,:), 'fontsize', fontsi);
    set(gca, 'xticklabel', '', 'yticklabel', '');
    if (m == 1)
        title ('Difference', 'fontsize', fontsi);
    end
    
    cbx = lmar+3*pwid+7*lrspc+clrspc;

    caxis = colorbar('location','manual','position',[cbx cby cwid chei]);
    caxiz(m,2) = caxis;
    set(caxis, 'fontsize', fontsi);
    if (m == 3)
        set(caxis, 'fontsize', fontsi);
        x1 = caxis.Position(1);
        y1 = caxis.Position(2);
        annotation(gcf, 'textbox', [x1+dxcl y1+dycl 0.2 0.07], 'string',{'Diff','(%)'}, 'linestyle', 'none', 'fontsize', fontsi);
        %xlabel(caxis, {'Diff','(%)'}, 'fontsize', fontsi);
    end
end

thesis_cmap;

for m = 1:3
    set(caxiz(m,2), 'ytick', [0 1 2]);
end

annotation(gcf,'textbox',[0.48 0.89 0.2 0.06],'Interpreter','latex','String',{'$\beta = 0$'}, 'fontsize', fontsi,'LineStyle','none','FontAngle','italic');
annotation(gcf,'textbox',[0.48 0.35 0.2 0.06],'Interpreter','latex','String',{'$\beta = \frac{\pi}{4}$'}, 'fontsize', fontsi,'LineStyle','none','FontAngle','italic');

for n = 1:length(axiz)
    multi_drawBodiesOnPlot(hbs, 'Color', 'k', 'Axis', axiz(n))
end

print('-dpng', '-r600', [folder 'pics\cyl_med']);

%% 5.31) Medium Array - Cylinder - check motions

%load cyl_med_array

Xiw = wcomp.Motions;
Xia = acomps.Motions;

pos = acomps.BodXY;

hind = 3:6:96;

Xiwh_b1 = squeeze(Xiw(2,1,hind));
Xiwh_b2 = squeeze(Xiw(2,2,hind));
Xiah_b1 = squeeze(Xia(2,1,hind));
Xiah_b2 = squeeze(Xia(2,2,hind));

figure;
subplot(2,2,1);
scatter3(pos(:,1), pos(:,2), abs(Xiwh_b1));
hold on;
scatter3(pos(:,1), pos(:,2), abs(Xiah_b1), 'Marker','.');

subplot(2,2,2);
scatter3(pos(:,1), pos(:,2), abs(Xiah_b1 - Xiwh_b1)./abs(Xiwh_b1));

subplot(2,2,3);
scatter3(pos(:,1), pos(:,2), abs(Xiwh_b2));
hold on;
scatter3(pos(:,1), pos(:,2), abs(Xiah_b2), 'Marker','.');

subplot(2,2,4);
scatter3(pos(:,1), pos(:,2), abs(Xiah_b2 - Xiwh_b2)./abs(Xiwh_b2));

%% 5.4) Medium Array - Attenuator - WAMIT

lam = [3, 10, 30];
beta = [0, pi/4];

rowSpc = 20;
lrSpc = 20;

pos = zeros(11,2);
%front row
pos(1:6,1) = -rowSpc/2*ones(6,1);
pos(1:6,2) = -2.5*lrSpc:lrSpc:2.5*lrSpc;
%back row
pos(7:11,1) = rowSpc/2*ones(5,1);
pos(7:11,2) = -2*lrSpc:lrSpc:2*lrSpc;

%pos = [0 -10; 0 10]

fieldArray = BemFieldArray([-30 -60 0], [0.5 0.5 1], [241 241 1]);

[wcomp, ~, ~, wWF] = multi_makeRun('Atten', 'MedArray', [folder 'runs'], lam, beta, pos, fieldArray);

wcompC = wcomp;

save([folder '\multi\atten_medArray'], 'wcomp', 'wWF');

% Angles - rotate front row
angles = zeros(11,1);
angles(1:6) = 45*ones(6,1);

[wcomp, ~, ~, wWF] = multi_makeRun('Atten', 'MedArray', [folder 'runs'], lam, beta, pos, fieldArray, 'Rotate', angles);

save([folder '\multi\atten_medArray_r'], 'wcomp', 'wWF');

%% 5.5) Medium Array - Attenuator - IT

% straight
load atten_medArray
wcomps = wcomp;
wWFs = wWF;
clear wcomp wWF
[X, Y] = wWFs.FieldPoints;

Nbod = length(wcomps.Bodies);

d = 35000;
Ndof = 8*Nbod;
D = zeros(Ndof,Ndof);
for n = 1:Nbod
    i7 = (n - 1)*8 + 7;
    i8 = i7 + 1;
    D(i7, i7) = d;
    D(i8, i8) = d;
end

wcomps.SetDpto(D);
wWFs.BodyMotions = wcomps.Motions;

load 'atten_hydBod';

for n = 1:Nbod
    hbs(n) = HydroBody(attenHydBod);
    hbs(n).XYpos = wcomps.Bodies(n).XYpos;
end

acomps = HydroArrayComp(hbs, wcomps.IncWaves);
acomps.SetDpto(D);

tic 
aWFs = acomps.WaveField(true, X, Y, 'NoVel');
toc

aWFs.BodyMotions = acomps.Motions;

% rotated
load atten_medArray_r
wcompr = wcomp;
wWFr = wWF;
clear wcomp wWF

wcompr.SetDpto(D);
wWFr.BodyMotions = wcompr.Motions;

for n = 1:Nbod
    hbs(n) = HydroBody(attenHydBod);
    hbs(n).XYpos = wcompr.Bodies(n).XYpos;
    hbs(n).Angle = wcompr.Bodies(n).Angle;
end

acompr = HydroArrayComp(hbs, wcompr.IncWaves);
acompr.SetDpto(D);

tic 
aWFr = acompr.WaveField(true, X, Y, 'NoVel');
toc

aWFr.BodyMotions = acompr.Motions;

save([folder '\multi\atten_medArray_WIT'], 'acomps', 'acompr', 'wcomps', 'wcompr', 'wWFs', 'wWFr', 'aWFs', 'aWFr')

%% 5.6) Medium Array - Attenuator - Plot

%load atten_medArray_WIT

hbs = acomps.Bodies;
[X, Y] = wWFs.FieldPoints;

fighei = 14;

type = 'Total';

clims = [0.7 1.3; 0.7 1.3; 0.7 1.3];
clims2 = [0 2; 0 2; 0 2];

etaW = wWFs.Elevation(type);
etaA = aWFs.Elevation(type);

etaWR = wWFr.Elevation(type);
etaAR = aWFr.Elevation(type);

figure
set(gcf, 'PaperPosition', [0 0 figwid fighei]);


pwid = 0.23;
phei = pwid;
lmar = 0.12;
lrspc = 0.013;
lrspc2 = 8.5*lrspc;
tbspc1 = 0.02;
tmar = 0.11;

chei = phei - 0.04;
cwid = 0.02;
clrspc = 0.02;
cbotspc = 0.02;

nb = 1;

axiz = zeros(9,2,3);
caxiz = zeros(3,2);
iax = 1;

for m = 1:3 
    
    nt = m;
    bpos = 1-tmar-(m-1)*(phei+tbspc1)-phei;
        
    if (m == 3)
        nb = 2;
        nt = 2;
        bpos = 1-tmar-(m-1)*(phei+2.5*tbspc1)-phei;
    end
        
    subplot('Position', [lmar bpos pwid phei]);
    if (m == 3)
        pcolor(X, Y, abs(etaWR{nt, nb}));
        multi_drawBodiesOnPlot(acompr.Bodies)
        axiz(iax, 3) = true; 
    else
        pcolor(X, Y, abs(etaW{nt, nb}));
        multi_drawBodiesOnPlot(acomps.Bodies)
        axiz(iax, 3) = false; 
    end
    axiz(iax, 1) = gca; 
    axiz(iax, 2) = false; 
    iax = iax + 1;
    fet;
    shading interp
    set(gca, 'clim', clims(m,:), 'fontsize', fontsi);
    if (m == 1)
        title ('WAMIT', 'fontsize', fontsi);
        ylabel('\lambda/a = 3', 'fontsize', fontsi);
        set(gca, 'xticklabel', '', 'yticklabel', '');
    elseif (m == 2)
        ylabel('\lambda/a = 10', 'fontsize', fontsi);
        set(gca, 'xticklabel', '', 'yticklabel', '');
    else
        xlabel('x/a', 'fontsize', fontsi);
        ylabel({'\lambda/a = 10','y/a'}, 'fontsize', fontsi);
        set(gca,'xtick', [-20 20 60], 'ytick', [-60 -20 20 60]);
    end

        
    subplot('Position', [lmar+pwid+lrspc bpos pwid phei]);
    if (m == 3)
        pcolor(X, Y, abs(etaAR{nt, nb}));
        multi_drawBodiesOnPlot(acompr.Bodies, 'Cir')
        axiz(iax, 3) = true; 
    else
        pcolor(X, Y, abs(etaA{nt, nb}));
        multi_drawBodiesOnPlot(acomps.Bodies, 'Cir')
        axiz(iax, 3) = false; 
    end
    axiz(iax, 1) = gca; 
    axiz(iax, 2) = true; 

    iax = iax + 1;
    fet;
    shading interp
    set(gca, 'clim', clims(m,:), 'fontsize', fontsi);
    set(gca, 'xticklabel', '', 'yticklabel', '');
    if (m == 1)
        title ('IT', 'fontsize', fontsi);
    end
    
    cbx = lmar+2*pwid+lrspc+clrspc;
    cby = bpos + cbotspc;

    caxis = colorbar('location','manual','position',[cbx cby cwid chei]);
    caxiz(m,1) = caxis;
    set(caxis, 'fontsize', fontsi);
    
    if (m == 3)        
        dxcl = -0.02;
        dycl = -0.08;
        x1 = caxis.Position(1);
        y1 = caxis.Position(2);
        annotation(gcf, 'textbox', [x1+dxcl y1+dycl 0.2 0.07], 'string', '|\eta/a|', 'linestyle', 'none', 'fontsize', fontsi);
    end
    
    subplot('Position', [lmar+2*pwid+lrspc2 bpos pwid phei]);
    %pcolor(X, Y, abs((etaA{m,nb}-etaW{m,nb})./etaW{m,nb}));
    if (m == 3)
        pcolor(X, Y, 100*abs((etaAR{nt,nb}-etaWR{nt,nb})));
        multi_drawBodiesOnPlot(acompr.Bodies, 'Cir')
        axiz(iax, 3) = true; 
    else
        pcolor(X, Y, 100*abs((etaA{nt,nb}-etaW{nt,nb})));
        multi_drawBodiesOnPlot(acomps.Bodies, 'Cir')
        axiz(iax, 3) = false; 
    end
    axiz(iax, 1) = gca; 
    axiz(iax, 2) = true; 
    iax = iax + 1;
    fet;
    shading interp
    set(gca, 'clim', clims2(m,:), 'fontsize', fontsi);
    set(gca, 'xticklabel', '', 'yticklabel', '');
    if (m == 1)
        title ('Difference', 'fontsize', fontsi);
    end
    
    cbx = lmar+3*pwid+lrspc2+clrspc;

    caxis = colorbar('location','manual','position',[cbx cby cwid chei]);
    caxiz(m,2) = caxis;
    set(caxis, 'fontsize', fontsi);
    if (m == 3)
        dxcl = -0.02;
        dycl = -0.08;
        x1 = caxis.Position(1);
        y1 = caxis.Position(2);
        annotation(gcf, 'textbox', [x1+dxcl y1+dycl 0.2 0.07], 'string', {'Diff', '(%)'}, 'linestyle', 'none', 'fontsize', fontsi);
    end
end

annotation(gcf,'textbox',[0.48 1-tmar+0.02 0.2 0.06],'Interpreter','latex','String',{'$\beta = 0$'}, 'fontsize', fontsi,'LineStyle','none','FontAngle','italic');
annotation(gcf,'textbox',[0.48 1-tmar-0.56 0.2 0.06],'Interpreter','latex','String',{'$\beta = \frac{\pi}{4}$'}, 'fontsize', fontsi,'LineStyle','none','FontAngle','italic');

thesis_cmap

print('-dpng', '-r600', [folder 'pics\atten_med']);

%% 6.1) Spectral Array - compute hydrobody

lam = 1:24;

[~, ~, attenHydBod] = multi_makeRun('Atten', 'HydroBody', [folder 'runs'], lam, []);

save([folder '\multi\atten_hydBod_spec'], 'attenHydBod');

%% 6.2) Spectral Array - compute wave field, power

load atten_hydBod_spec

% make spectrum
Hs = 1;
lamp = 10;
h = attenHydBod.H;

Tp = IWaves.Lam2T(lamp,h);
fp = 1./Tp;

T = attenHydBod.T;
f = 1./T;

s = 10;
thetac = 0;

beta = [-pi/4 -pi/8 0 pi/8 pi/4];

Spec = Bretschneider(Hs, Tp, T, 'Cos2s', s, thetac, beta);

% make array from medium array
load atten_medArray
[X, Y] = wWF.FieldPoints;

Nbod = length(wcomp.Bodies);

d = 35000;
Ndof = 8*Nbod;
D = zeros(Ndof,Ndof);
for n = 1:Nbod
    i7 = (n - 1)*8 + 7;
    i8 = i7 + 1;
    D(i7, i7) = d;
    D(i8, i8) = d;
    
    hbs(n) = HydroBody(attenHydBod);
    hbs(n).XYpos = wcomp.Bodies(n).XYpos;
end

acomp = HydroArrayComp(hbs, Spec.GetPlaneWaves(h));
acomp.SetDpto(D);

tic 
aWF = acomp.WaveField(true, X, Y, 'NoVel');
toc

save([folder '\multi\atten_array_spec'], 'acomp', 'aWF', 'Spec', '-v7.3');

%% 6.3) Spectral Array - compute q factor, Hs/Hs0

load atten_array_spec

rho = 1000;
h = 10;

powa = acomp.Power;
Nbod = length(acomp.Bodies);
powA = zeros(Nbod,1);

for n = 1:Nbod
    i7 = (n - 1)*8 + 7;
    i8 = i7 + 1;
    
    p7 = squeeze(powa(:,:,i7));
    p8 = squeeze(powa(:,:,i8));
    
    powA(n) = sum(sum(p7)) + sum(sum(p8));
end

load atten_hydBod_spec;

bcomp = HydroBodyComp(attenHydBod, Spec.GetPlaneWaves(h));
d = 35000;
Ndof = 8;
D = zeros(Ndof,Ndof);
D(7,7) = d;
D(8,8) = d;
bcomp.SetDpto(D);

powb = bcomp.Power;

p7 = squeeze(powb(:,:,7));
p8 = squeeze(powb(:,:,8));

powB = sum(sum(p7)) + sum(sum(p8));

pos = acomp.BodXY;
qA = powA./powB;

Hs = aWF.SigWaveHeight('Total');
Hs0 = Spec.SigWaveHeight;

[X, Y] = aWF.FieldPoints;

%% 6.4) Spectral Array - plot

fighei = 16;

figure
set(gcf, 'PaperPosition', [0 0 figwid fighei]);


pwid1 = 0.6;
pwid2 = 0.8;

phei1 = 0.4;
phei2 = 0.4;

lmar1 = 0.22;
lmar2 = 0.07;
lrspc = 0.12;
tbspc = 0.02;
tmar = 0.05;

% chei = phei - 0.04;
% cwid = 0.02;
% clrspc = 0.02;
% cbotspc = 0.02;

indsp1 = find(qA >= 1);
c1 = [0.9290, 0.6940, 0.1250];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0, 0.4470, 0.7410];
col = ones(length(qA),1)*c2;
col(indsp1,:) = ones(length(indsp1),1)*c3;

siz = 50;

axq = subplot('Position', [lmar1 1-tmar-phei1+0.05 pwid1 phei1]);

scatter3(pos(:,1), pos(:,2), qA, siz, col, 'Marker','.');
x10s = [-30 30 30 -30];
y10s = [-60 -60 60 60];
[X10s, Y10s] = meshgrid(x10s, y10s);
bod1 = ones(size(X10s));

hold on;
pla = surf(X10s, Y10s, bod1);
alpha(pla, 0.1);

%view(gca,[-6 2]);
view(gca,[-7 5]);
set(gca, 'xlim', [-30 30], 'ylim', [-60 60], 'zlim', [0 1.2], 'xtick', [-20 -10 0 10 20 30], 'ztick', [0 0.5 1 1.2], 'zticklabel', {'0', '0.5', '1', ''}, 'fontsize', fontsi);
set(gca, 'dataaspectratio', [1 1 1/20])
set(gca, 'plotboxaspectratio', [1 1 1])

xlabel('x/d', 'fontsize', fontsi);
ylabel('y/d', 'fontsize', fontsi);
zlabel('q    ', 'fontsize', fontsi);
set(get(gca, 'zlabel' ), 'Rotation' ,0 );

subplot('Position', [lmar2 1-tmar-phei1-tbspc-phei2 pwid2 phei2]);

pcolor(X, Y, Hs./Hs0);
fet;
shading interp;
thesis_cmap;
set(gca, 'clim', [0.9 1.1])
multi_drawBodiesOnPlot(acomp.Bodies, 'Cir')

set(gca,'xtick', [-20 0 20 40 60 80], 'ytick', [-50 -30 -10 10 30 50], 'fontsize', fontsi);
xlabel('x/d', 'fontsize', fontsi);
ylabel('y/d', 'fontsize', fontsi);

axpos = get(gca, 'position');

caxis = colorbar;

chei = axpos(4) - 0.04;
cwid = 0.05;
clrspc = -0.08;
cbotspc = 0.02;

cbx = axpos(1) + axpos(3) + clrspc;
cby = axpos(2) + cbotspc;

%caxis = colorbar('location','manual','position',[cbx cby cwid chei]);
set(caxis, 'Position', [cbx cby cwid chei], 'fontsize', fontsi);
dxcl = -0.02;
dycl = -0.08;
x1 = caxis.Position(1);
y1 = caxis.Position(2);
annotation(gcf, 'textbox', [x1+dxcl y1+dycl 0.2 0.07], 'string', 'H_{s}/H^{I}_{s}', 'linestyle', 'none', 'fontsize', fontsi);

xst = 0.45;
yst1 = 0.55;
yst2 = 0.0;
annotation(gcf,'textbox',[xst yst1 0.3 0.05],'String',{'a) Interaction factor'}, 'fontsize', fontsi,'LineStyle','none');
annotation(gcf,'textbox',[xst yst2 0.3 0.05],'String',{'b) Wave field'}, 'fontsize', fontsi,'LineStyle','none');


print('-dpng', '-r600', [folder 'pics\spec_array']);

%% 7.1) Large Array - hydrobody 
% recompute a hydrobody at just lam = 10;

lam = 10;

[~, ~, attenHydBod] = multi_makeRun('Atten', 'HydroBody', [folder 'runs'], lam, []);

save([folder '\multi\atten_hydBod_large'], 'attenHydBod');

%% 7.2) Large Array - compute

%load atten_hydBod_large

Nbod = 101;

rowSpc = 20;
lrSpc = 20;

pos = zeros(Nbod,2);
row1 = 34;
row2 = 33;
row3 = 34;
%first row
pos(1:row1,1) = -rowSpc*ones(row1,1);
pos(1:row1,2) = -(row1-1)/2*lrSpc:lrSpc:(row1-1)/2*lrSpc;
%second row
pos(row1+1:row1+row2,1) = zeros(row2,1);
pos(row1+1:row1+row2,2) = -(row2-1)/2*lrSpc:lrSpc:(row2-1)/2*lrSpc;
%third row
pos(row1+row2+1:row1+row2+row3,1) = rowSpc*ones(row3,1);
pos(row1+row2+1:row1+row2+row3,2) = -(row3-1)/2*lrSpc:lrSpc:(row3-1)/2*lrSpc;

% q-factor
beta = 0;

load 'atten_hydBod_large';

attenHydBodL = attenHydBod;

for n = 1:Nbod
    hbs(n) = HydroBody(attenHydBodL);
    hbs(n).XYpos = pos(n,:);
end

acomp = HydroArrayComp(hbs, PlaneWaves(1, attenHydBodL.T, beta, attenHydBodL.H));

d = 35000;
Ndof = 8*Nbod;
D = zeros(Ndof,Ndof);
for n = 1:Nbod
    i7 = (n - 1)*8 + 7;
    i8 = i7 + 1;
    D(i7, i7) = d;
    D(i8, i8) = d;
end

acomp.SetDpto(D);

pA = squeeze(acomp.Power);

powA = zeros(Nbod,1);
for n = 1:Nbod
    i7 = (n - 1)*8 + 7;
    i8 = i7 + 1;
    powA(n) = pA(i7) + pA(i8);
end

load atten_hydBod;

bcomp = HydroBodyComp(attenHydBod, PlaneWaves(1, attenHydBod.T, beta, attenHydBod.H));
d = 35000;
Ndof = 8;
D = zeros(Ndof,Ndof);
D(7,7) = d;
D(8,8) = d;
bcomp.SetDpto(D);

powb = bcomp.Power;
powB = powb(2,1,7) + powb(2,1,8);

pos = acomp.BodXY;
qA = powA./powB;

% wave field
x = -40:1.4:100;
y = -350:1.4:350;

[X, Y] = meshgrid(x, y);

tic 
aWF = acomp.WaveField(true, X, Y, 'NoVel');
wftime = toc

etaA = aWF.Elevation('Total');

save([folder '\multi\atten_array_large'], 'acomp', 'qA', 'pos', 'X', 'Y', 'aWF', 'etaA', '-v7.3');

%% 7.3) Large Array - plot

fighei = 22;

figure
set(gcf, 'PaperPosition', [0 0 figwid fighei]);


pwid1 = 0.9;
pwid2 = 0.9;

phei1 = 0.2;
phei2 = 0.55;

lmar = 0.11;
lrspc = 0.1;
tbspc = 0.1;
tmar = 0.04;

% chei = phei - 0.04;
% cwid = 0.02;
% clrspc = 0.02;
% cbotspc = 0.02;

indsp1 = find(qA >= 1);
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0, 0.4470, 0.7410];
col = ones(length(qA),1)*c2;
col(indsp1,:) = ones(length(indsp1),1)*c3;

% col = zeros(size(qA));
% col(indsp1) = ones(size(indsp1));

siz = 50;

axq = subplot('Position', [lmar 1-tmar-phei1 pwid1 phei1]);
scatter3(pos(:,1), pos(:,2), qA, siz, col, 'Marker','.');
x10s = [-40 40 40 -40];
y10s = [-350 -350 350 350];
[X10s, Y10s] = meshgrid(x10s, y10s);
bod1 = ones(size(X10s));

hold on;
pla = surf(X10s, Y10s, bod1);
alpha(pla, 0.1);

view(gca,[-1.5 1]);
%view(gca,[-7 5]);
set(gca, 'xlim', [-40 40], 'ylim', [-350 350], 'zlim', [0 1.5], 'xtick', [-20 0 20 40], 'ytick', [-200 0 200], 'fontsize', fontsi);
set(gca, 'dataaspectratio', [1 1 1/30])
set(gca, 'plotboxaspectratio', [1 1 1])

%title('\lambda/a = 10, \beta = 0','fontname',font, 'fontsize', fontsi, 'fontweight', fontwe);
xlabel('x/d', 'fontsize', fontsi);
ylabel('y/d   ', 'fontsize', fontsi);
zlabel('q', 'fontsize', fontsi);
set(get(gca, 'zlabel' ), 'Rotation' ,0 );

subplot('Position', [lmar 1-tmar-phei1-tbspc-phei2 pwid2 phei2]);

pcolor(X, Y, abs(etaA{1}));
fet;
shading interp;
thesis_cmap;
set(gca, 'clim', [0.7 1.3])
multi_drawBodiesOnPlot(acomp.Bodies)

set(gca, 'fontsize', fontsi);
xlabel('x/d', 'fontsize', fontsi);
ylabel('y/d', 'fontsize', fontsi);

caxis = colorbar;
set(caxis, 'fontsize', fontsi);
dxcl = 0.02;
dycl = -0.08;
x1 = caxis.Position(1);
y1 = caxis.Position(2);
annotation(gcf, 'textbox', [x1+dxcl y1+dycl 0.2 0.07], 'string', '|\eta/a|', 'linestyle', 'none', 'fontsize', fontsi);
%xlabel(caxis, '|\eta/a|', 'fontsize', fontsi);

xst = 0.4;
yst = 0.01;
annotation(gcf,'textbox',[xst 0.65 0.3 0.06],'String',{'a) Interaction factor'}, 'fontsize', fontsi,'LineStyle','none');

annotation(gcf,'textbox',[xst 0 0.3 0.06],'String',{'b) Wave field'}, 'fontsize', fontsi,'LineStyle','none');

print('-dpng', '-r600', [folder 'pics\large_array']);
