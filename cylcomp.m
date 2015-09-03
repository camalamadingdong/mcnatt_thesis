%% Cylindrical Surface Computation
% copied from ewtec_PaperPics in N:\RDS\WAMIT\wf_study\ewtec_2013

folder = 'C:\Users\s1213969\Dropbox\matlab\thesis\';

fontsi = 10;
figwid = 14;

d = 1;
h = 10*d;
rho = 1000;
g = IWaves.G;

%% Plot the geometries


fighei = 5;

dwid = 0.12;
wids = [0.1 0.2 0.38];
hei = 0.7;

lmar = 0.08;
lefs = lmar*[1 1 1] + dwid*[0 1 2] + wids(1)*[0 1 1] + wids(2)*[0 0 1] + [0 0 -0.05];
bot = 0.1;

ytit = 0.8;

figure;
set(gcf, 'PaperPosition', [0 0 figwid fighei]);

cyl = cylcomp_createGeo(d, 'Cyl');
subplot('Position', [lefs(1), bot, wids(1), hei]);
%plot(cyl.PanelGeo);
plot(cyl.PanelGeo, 'OnlyWet');
set(gca, 'view', [-37.5000 30]);
axis equal; 
set(gca, 'xtick', [], 'ytick', [], 'ztick', [-1.5 0], 'FontSize', fontsi);
% xlabel('d', 'FontSize', fontsi);
% ylabel('d', 'FontSize', fontsi);
% zlabel('c_c', 'FontSize', fontsi);
xtit = lefs(1)-0.02;
annotation(gcf, 'textbox', [xtit ytit 0.1 0.1], 'string', 'Cylinder', 'linestyle', 'none', 'fontweight', 'bold', 'fontsize', fontsi);
%title(['Cylinder',''], 'FontSize', fontsi);


flap = cylcomp_createGeo(d, 'Flap');
subplot('Position', [lefs(2), bot, wids(2), hei]);
plot(flap.PanelGeo, 'OnlyWet');
axis equal; 
set(gca, 'xtick', [], 'ytick', [-2 0 2], 'ztick', [-1 0], 'FontSize', fontsi);
xtit = lefs(2)+0.07;
annotation(gcf, 'textbox', [xtit ytit 0.1 0.1], 'string', 'Flap', 'linestyle', 'none', 'fontweight', 'bold', 'fontsize', fontsi);
% xlabel('d_f', 'FontSize', fontsi);
% ylabel('a_f', 'FontSize', fontsi);
% zlabel('c_f', 'FontSize', fontsi);
%title('Flap', 'FontSize', fontsi);


atten = cylcomp_createGeo(d, 'Atten');
subplot('Position', [lefs(3), bot, wids(3), hei]);
plot(atten.PanelGeo, 'OnlyWet');
axis equal;
set(gca, 'xtick', [-4 -2 0 2 4], 'ytick', [], 'ztick', [], 'FontSize', fontsi);
xtit = lefs(3)+0.12;
annotation(gcf, 'textbox', [xtit ytit 0.1 0.1], 'string', 'Attenuator', 'linestyle', 'none', 'fontweight', 'bold', 'fontsize', fontsi);
% xlabel('l_a', 'FontSize', fontsi);
% ylabel('d_a', 'FontSize', fontsi);
%title(['Attenuator',''], 'FontSize', fontsi);

print('-dpng', '-r400', [folder 'pics\cylcomp_geos']);

%% Plot WEC with surf points


fighei = 12;

atten = cylcomp_createGeo(d, 'Atten');
figure;
set(gcf, 'PaperPosition', [0 0 figwid fighei]);
plot(atten.PanelGeo, 'OnlyWet');
axis equal
axis off;
%set(gca, 'FontSize', fontsi, 'ytick', [-0.2 0 0.2], 'yTickLabel',{'', '', '0.2'});

dr = 0.1;
r = 5*d;
r = r + dr*d;

nZ = 21;
nTheta = 2^5;

z = -h*(1-cos(0:pi/2/nZ:pi/2));
z = round(z*10^6)/10^6;
points = makeCirWFPoints(r, nTheta, z);

hold on

s = 50;
co = [0 0.5 0];
scatter3(points(:,1), points(:,2), points(:,3), s, co, 'Marker','.');

print('-dpng', '-r400', [folder 'pics\cylcomp_surfPoints']);

%% Performance

figwid = 10;
fighei = 10;

beta = 0;
lam = 1:0.5:50;

wftype = 'None';

[comph] = cylcomp_comp([folder '\runs'], d, lam, beta, 'Heave', wftype);
[comps] = cylcomp_comp([folder '\runs'], d, lam, beta, 'Surge', wftype);
[compf] = cylcomp_comp([folder '\runs'], d, lam, beta, 'Flap', wftype);
[compa] = cylcomp_comp([folder '\runs'], d, lam, beta, 'Atten', wftype);

Ph = comph.Power;
Ps = comps.Power;
Pf = compf.Power;
Pa = compa.Power;
Pa = Pa(:,2);

T = IWaves.Lam2T(lam,h);
Ef = IWaves.UnitEnergyFlux(rho, T, h).';

Ph = Ph./Ef;
Ps = Ps./Ef;
Pf = Pf./Ef;
Pa = Pa./Ef;

figure;
set(gcf, 'PaperPosition', [0 0 figwid fighei]);

plot(lam, [Ph, Ps, Pf, Pa]);

set(gca, 'xlim', [0 40], 'ylim', [0 1.2], 'fontsize', fontsi);
xlabel('Wavelength (\lambda/d)','fontsize', fontsi);
ylabel('Capture Width (CW/d)','fontsize', fontsi);
legend('Heave', 'Surge', 'Flap', 'Attenuator');

print('-dpng', '-r400', [folder 'pics\cylcomp_perf']);

%% M values plot 

lam = 10;
beta = 0;

type = 'Heave';
wfdist = 20;

wftype = 'Grid';
[~, waveG] = cylcomp_comp([folder '\runs'], d, lam, beta, type, wftype, wfdist);

wftype = 'Circ';
[comp, waveSurf] = cylcomp_comp([folder '\runs'], d, lam, beta, type, wftype, wfdist);

M = 8;
N = 20; 

[X Y] = waveG.FieldPoints;

ErrPlot = cell(1,M+1);

waveG.BodyMotions = comp.Motions;
eta = waveG.Elevation('Total');
etaG = eta{1};

for m = 0:M

    [waveC2, waveF, r0] = cylcomp_makeCirWF(waveSurf, m, N, X, Y);
    waveC2.BodyMotions = comp.Motions;
    
    eta = waveC2.Elevation('Total');
    etaC = eta{1};
    
    err = abs(etaC - etaG)./abs(etaG);

    ErrPlot{1+m} = err;
end

scale = 100;

fighei = 16;

tmar = 0.07;
lmar = 0.12;
lrspc = 0.05;
tbspc = 0.05;

pwid = 0.16;
phei = 0.16;

figure;
set(gcf, 'PaperPosition', [0 0 figwid fighei]);

clim = scale*[0 0.01];

R = sqrt(X.^2 + Y.^2);
inds = (R > r0);
Err = zeros(1,M+1);
Npts = sum(sum(inds));

thetar = 0:2*pi/100:2*pi;
cirx = r0*cos(thetar);
ciry = r0*sin(thetar);

Pleft = [lmar, lmar+pwid+lrspc, lmar+2*pwid+2*lrspc];
Pbot = [1-tmar-phei, 1-tmar-2*phei-tbspc, 1-tmar-3*phei-2*tbspc];

for m = 0:M
    
    err = ErrPlot{1+m};
    Err(1+m) = sum(err(inds))./Npts;
    
    pleft = Pleft(rem(m,3)+1);
    pbot = Pbot(floor(m/3)+1);
    
    subplot('position', [pleft pbot pwid phei])
    
    pcolor(X,Y,scale*err);
    
    shading interp;
    axis equal;
    axis tight;
    set(gca, 'clim', clim)
    title(['N = ' num2str(m)], 'fontsize', fontsi);
    
    set(gca, 'fontsize', fontsi, 'xlim', [-10 10], 'ylim', [-10 10]);
    if (m == 6)
        xlabel('x/d', 'fontsize', fontsi);
        ylabel('y/d', 'fontsize', fontsi);
    else
        set(gca, 'xtick', [], 'ytick', []);
    end
    
    hold on;
    plot(cirx, ciry, 'w');
end

cmap = colormap(gca);
cmap2 = zeros(10,3);
for n = 1:10
    cmap2(n,:) = cmap(1+(n-1)*7,:); 
end
colormap(cmap2)

cbx = 0.78;
cby = 0.385;
chei = 0.5;
cwid = 0.05;

caxis = colorbar('position',[cbx cby cwid chei]);
%caxis = colorbar('location','manual','position',[cbx cby cwid chei],'plotboxaspectratio',[1 10 1]);
set(caxis, 'fontsize', fontsi);
xlabel(caxis, 'Err (%)', 'fontsize', fontsi);

subplot('position', [0.12 0.1 0.8 0.175])

plot(0:M,scale*Err,'MarkerFaceColor',[0 0 1], 'MarkerEdgeColor',[0 0 1],'Marker','o', 'LineStyle','none','MarkerSize',2);
set(gca, 'fontsize', fontsi);
xlabel('N','fontsize', fontsi);
ylabel('Mean Err (%)', 'fontsize', fontsi);


print('-dpng', '-r400', [folder 'pics\cylcomp_Ms']);

%% Make and save all wave fields

M = 20;
wfdist = 20;
wftype = 'Both';

lams = [3 10 30];
betas = [0 pi/6];
bstrs = {'0', '30'};
types = {'Heave', 'Surge', 'Flap', 'Atten2'};

tic
for l = 2:3
    lam = lams(l);
    
    for m = 1:2
        beta = betas(m);
        bstr = bstrs{m};
        
        for n = 4:4
            type = types{n};
            
            [comp, waveG, waveSurf] = cylcomp_comp([folder '\runs'], d, lam, beta, type, wftype, wfdist);

            N = cylcomp_checkEvCut(waveSurf, waveG, comp);

            [X, Y] = waveG.FieldPoints;
            [waveC, waveF, r0, B] = cylcomp_makeCirWF(waveSurf, N, M, X, Y);
            
            nameStr = ['cylc_' type '_lam' num2str(lam) '_b' bstr];
            
            save([folder 'cylcomp\' nameStr], 'nameStr', 'waveG', 'waveC', 'waveF', 'comp', 'r0', 'B', 'M', 'N', 'wfdist', 'lam', 'beta', 'type');
            
            if (l == 3)
                clim1 = [0.98 1.02];
                clim2 = [0 1];
                clim3 = [0 10];
                
                cylcomp_plotWF([folder 'pics\' nameStr], waveG, waveC, waveF, comp, r0, B, 'clims', clim1, clim2, clim3);
            else
                cylcomp_plotWF([folder 'pics\' nameStr], waveG, waveC, waveF, comp, r0, B);
            end
            
            close all
        end
    end
end
toc

%% Plot WF test

lam = 10;
beta = 0;
wfdist = 20;
M = 20;


type = 'Flap';
wftype = 'Circ';
[~, waveSurf] = cylcomp_comp([folder '\runs'], d, lam, beta, type, wftype, wfdist);

wftype = 'Grid';
[comp, waveG] = cylcomp_comp([folder '\runs'], d, lam, beta, type, wftype, wfdist);

N = cylcomp_checkEvCut(waveSurf, waveG, comp);

[X, Y] = waveG.FieldPoints;
[waveC, waveF, r0, B] = cylcomp_makeCirWF(waveSurf, N, M, X, Y);


cylcomp_plotWF([folder 'pics\test1'], waveG, waveC, waveF, comp, r0, B)

%% Find far-field error distance

lam = 3;
dir = 0;
type = '4_a'; 

name = ['l' num2str(lam) 'b' num2str(dir) '_' type];
load(name);

eta = waveG.Elevation('Total');
etaG = eta{1};

eta = waveF.Elevation('Total');
etaF = eta{1};

Err = 100*abs(etaF-etaG)./abs(etaG);

[X Y] = waveG.FieldPoints;

R = sqrt(X.^2 + Y.^2);

endr = 100;
r0 = 0;
lim = 1;

for r = 0:1:endr
    inds = (R > r);
    err = Err(inds);
    if (all(err < lim))
        r0 = r;
        break;
    end
end

figure;
clim = [0 5];
pcolor(X,Y,Err);
shading interp;
axis equal;
axis tight;
set(gca, 'clim', clim)

thetar = 0:2*pi/100:2*pi;
cirx = r0*cos(thetar);
ciry = r0*sin(thetar);

hold on;
plot(cirx, ciry, 'w');

% colormap
cmap = colormap(gca);
cmap2 = zeros(10,3);
for n = 1:10
    cmap2(n,:) = cmap(1+(n-1)*7,:); 
end
colormap(cmap2)

display(['Error is less than ' num2str(lim) '% at r = ' num2str(r0/lam) '*lambda']);


