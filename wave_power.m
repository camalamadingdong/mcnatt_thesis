%% Wave power computation
% copied from ewtec_pap2 in Dropbox\matlab\studies
% 0) set up stuff

folder = 'C:\Users\s1213969\Dropbox\matlab\thesis\';
fontsi = 10;
figwid = 14;

%% 1) Bessel function pic

fighei = 6;

m1 = [0 1];
m2 = [7 12];
kr1 = linspace(0,15,101);

limM1 = sqrt(m1 + [1 1]);
limM2 = sqrt(m2 + [1 1]);

[M1, KR] = meshgrid(m1, kr1);
[M2] = meshgrid(m2, kr1);

J1 = besselj(M1, KR);
Y1 = bessely(M1, KR);
J2 = besselj(M2, KR);
Y2 = bessely(M2, KR);


figure;
set(gcf,'PaperPosition',[0 0 figwid fighei]);

dy = 2*10^-5;

subplot(1,2,1);
l0 = plot(kr1,J1);
hold on;
l2 = plot(kr1,Y1,'--');
set(l2(1), 'Color', l0(1).Color)
set(l2(2), 'Color', l0(2).Color)
%plot([1; 1]*limM1, [-2 -2; 2 2], 'k:')
set(gca, 'ylim', 2.2*[-1 1], 'ytick', [-2 -1 0 1 2],'FontSize', fontsi);
leg = legend('J_0', 'J_1', 'Y_0', 'Y_1');
%pos = leg.Position;
pos = [0.32 0.69 0.12 0.15];
set(leg, 'Box', 'off', 'Position', pos, 'FontSize', 6.5);
%ylabel('Magnitude');
xlabel('kr','FontSize', fontsi);

subplot(1,2,2);
l1 = plot(kr1,J2);
hold on;
l2 = plot(kr1,Y2,'--');
set(l2(1), 'Color', l0(1).Color)
set(l2(2), 'Color', l0(2).Color)
%plot([1; 1]*limM2, [-2 -2; 2 2], 'k:')
set(gca, 'ylim', 2.2*[-1 1], 'ytick', [-2 -1 0 1 2],'FontSize', fontsi);
leg = legend('J_7', 'J_{12}', 'Y_7', 'Y_{12}');
pos = [0.76 0.69 0.12 0.15];
set(leg, 'Box', 'off', 'Position', pos, 'FontSize', 6.5);
%ylabel('Magnitude');
xlabel('kr','FontSize', fontsi);

print('-dpng', '-r400', [folder 'pics\bessel_pow']);

%% 2) Point absorber pic
fighei = 7.5;

beta = pi/4;
M = 3;
m = -M:M;
a = 1;
aI = a/2*exp(-1i*m*(beta+pi/2));

close all;
figure;
set(gcf, 'PaperPosition', [0 0 figwid fighei]);
plotCoefs(M, aI, true, 'Color', [1 0 0], 'LineWidth', 1)

bR = a/2*[-1i*exp(1i*beta) -1 1i*exp(-1i*beta)];

MR = 1;

plotCoefs(MR, bR, false, 'Color', [0 0 1], 'LineWidth', 1, 'LineStyle', '--')

ax = gca;
lines = ax.Children;
leg = legend([lines(20), lines(1)], {'1/2 Incident', 'Radiated'}, 'Position', [0.7 0.7 0.1 0.1], 'FontSize', fontsi);
legend('boxoff');

x = 0.7;
y = 0.55;
w = 0.1;
h = 0.1;
an = annotation('textbox', [x y w h], 'String', ['\beta = \pi/4'], 'LineStyle', 'none', 'FontSize', fontsi);

% add labels
x = 0.52;
y = 0.55;
w = 0.1;
h = 0.1;
an = annotation('textbox', [x y w h], 'String', '|a|', 'LineStyle', 'none', 'FontSize', fontsi);
set(an, 'LineStyle', 'none');
x = 0.42;
y = 0.4;
an = annotation('textbox', [x y w h], 'String', 'ang(a)', 'LineStyle', 'none', 'FontSize', fontsi);
set(an, 'LineStyle', 'none');

% slope = [0.27/M -0.137/M];
slope = [0.26/M -0.16/M];
start = [0.5 0.38];

for m = -M:M
    x = start(1) + slope(1)*m;
    y = start(2) + slope(2)*m;
    if (m == 0)
        x = x + 0.02;
        y = y + 0.07;
        mstr = 'm = 0';
    else
        mstr =  num2str(m);
    end
    annotation('textbox', [x y w h], 'String', mstr, 'LineStyle', 'none', 'FontSize', fontsi);
end

print('-dpng', '-r400', [folder 'pics\rad_coefs']);

%% 3.1) Compute values for different beams for main pic

beams = [1 2 5 10 20 50];

Nb = length(beams);
M = 64;
nTheta = M*2^3;  

aSs = cell(Nb, 1);
aRs = cell(Nb, 2);

for n = 1:Nb
    [aS, aR] = wave_power_wamComp([folder 'runs'], beams(n), nTheta, M);
    aSs{n} = aS;
    aRs(n,1:2) = aR;
end

save([folder '\wave_power\pow_beams'], 'beams', 'aSs', 'aRs', 'M', 'Nb');

%% 3.2) Compute values for large range of beams pic

beamsB = 1:50;

NbB = length(beamsB);
M = 64;
nTheta = M*2^3;  

aSsB = cell(NbB, 1);
aRsB = cell(NbB, 2);

for n = 1:NbB
    [aS, aR] = wave_power_wamComp([folder 'runs'], beamsB(n), nTheta, M);
    aSsB{n} = aS;
    aRsB(n,1:2) = aR;
end

save([folder '\wave_power\pow_big_beams'], 'beamsB', 'aSsB', 'aRsB', 'M', 'NbB');

%% 4.1) Make main beams figure

load pow_beams

fighei = 18;
fontsi = 8;

M = 64;
Mt = 20;
inds = M+1-Mt:M+1+Mt;

m = -Mt:Mt;
aI = exp(-1i*(-M:M)*pi/2);
aI = aI.';

r = 1;
dep = 3;

h = 10;          
lam = 10;
T = IWaves.Lam2T(lam, h);
t = linspace(0,T,101);
omega = 2*pi/T;
k = 2*pi/lam;
eiomegat = exp(1i*omega*t);
wavx = linspace(-15,15,101);
wavz = exp(-1i*k*wavx);

rho = 1000;
g = 9.806650;
cg = IWaves.GroupVel(T, h);
const = 2*rho*g*cg/k;

uef = IWaves.UnitEnergyFlux(rho, T, h);
uef = uef./const;

figure;
set(gcf,'PaperPosition',[0 0 figwid fighei]);

tmar = 0.05;
lmar = 0.05;

fighei = 0.13;
%figwid1 = 0.14;
figwid1 = 0;
figwid2 = 0.2;
figwid3 = 0.25;

tbspc = 0.015;
lrspc = 0.055;

nsub = 5;

sca = [0.3 0.8 1.5 1.7 1.7 1.7];

bstartx = 0.01;
txtstarty = 1-tmar-0.072;

annotation('textbox', [bstartx txtstarty 0.1 0.1], 'String',{'Beam', '(B/A)'},...
        'LineStyle', 'none', 'FontWeight', 'bold', 'FontSize', fontsi)
% annotation('textbox', [0.11 txtstarty 0.1 0.1], 'String','Cylinder View',...
%         'LineStyle', 'none', 'FontWeight', 'bold', 'FontSize', fntsiz)
annotation('textbox', [0.12 txtstarty 0.1 0.1], 'String','Orbital Motion',...
        'LineStyle', 'none', 'FontWeight', 'bold', 'FontSize', fontsi)
cwstartx = lmar + 3*lrspc + figwid1 + figwid2 + 2*figwid3+0.01;
annotation('textbox', [cwstartx+0.01 txtstarty 0.1 0.1], 'String',{'E'},...
    'LineStyle', 'none', 'FontWeight', 'bold', 'FontSize', fontsi)


xs = [0.175 0.465 0.773 0.94];
vals = {'a)', 'b)', 'c)', 'd)'};
for n = 1:4
    annotation('textbox', [xs(n) 0.03 0.1 0.01], 'String', vals{n},...
        'LineStyle', 'none', 'FontWeight', 'bold', 'FontSize', fontsi)
end

for n = 1:Nb
    
    aS = aSs{n}.';
    aD = 1/2*aI + aS;
    aR1 = aRs{n,1}.';
    aR2 = aRs{n,2}.';
    BR = [aR1 aR2];
    zeta = BR\(-aD);
    d = aD + BR*zeta;
    
    pow = (abs(aD).^2 - abs(d).^2);
    powDiff = abs(aD).^2;
    eff = pow./powDiff;
    
    totPow = sum(pow);
    CW = totPow/uef;
    

    % Beam
    yplot = 1-tmar-fighei-(n-1)*(fighei+tbspc);
    beam = beams(n);
    annotation('textbox', [bstartx+0.01 yplot-0.01 0.1 0.1], 'String',...
        {[num2str(beam)]},...
        'LineStyle', 'none', 'FontSize', fontsi)
    
    % Plot motions
    xpath = real(zeta(1)*eiomegat);
    zpath = real(zeta(2)*eiomegat) - dep*ones(size(eiomegat));
    % arrow points
    xa = [xpath(1), xpath(25), xpath(50), xpath(75)];
    za = [zpath(1), zpath(25), zpath(50), zpath(75)];
    
    ua = omega*(za+dep*ones(size(za)));
    va = -omega*xa;
    la = sqrt(ua.^2 + va.^2);
    ua = ua./la;
    va = va./la;
    %subplot(Nb,nsub,nsub*(n-1)+2);
    subplot('Position', [lmar+lrspc+figwid1 yplot-0.0175 figwid2-0.018 fighei]);
    plot(xpath, zpath, 'Color', [0 0 0], 'LineWidth', 1);
    hold on;
    quiver(xa, za, ua, va, sca(n), 'k');
    plot(wavx, wavz, 'k--');
    %plot([-30 30], [0 0], 'k:');
    plot([-30 30], [-h -h], 'k');
    axis equal
    axis off
    set(gca, 'xlim', [-15 15], 'ylim', [-20 4]);

%     annotation('textbox', [lmar+lrspc+figwid+0.02 yplot+0.04 0.1 0.1], 'String',...
%         {['\zeta_1 = ' num2str(abs(zeta(1)), '%4.2f') ', \zeta_2 = ' num2str(abs(zeta(2)), '%4.2f')]},...
%         'LineStyle', 'none', 'FontSize', fntsiz)

    %subplot(Nb,nsub,nsub*(n-1)+4);
    subplot('Position', [lmar+3*lrspc+figwid1+figwid2+figwid3 yplot figwid3 fighei]);
    stem(m, eff(inds), 'Marker', '.', 'Color', [0 0 0]);
    hold on;
    plot([-Mt Mt], [1 1], 'k--');
    set(gca, 'ylim', [-0.1 1.1],'FontSize', fontsi);
    
    if (n == 1)
        title('Efficiency', 'FontSize', fontsi);
    end
    
    if (n == Nb)
        xlabel('m');
        set(gca, 'xtick', [-20 -10 0 10 20]);
    else
        set(gca, 'xtick', []);
    end
    CW2 = CW*k;
    annotation('textbox', [cwstartx+0.005 yplot-0.01 0.1 0.1], 'String',...
        {[num2str(CW2, '%3.1f')]}, 'LineStyle', 'none','FontSize', fontsi)
    plot([CW2/2 CW2/2], [0 1], 'k','LineStyle',':','Color', [1 0 0]);
    plot(-[CW2/2 CW2/2], [0 1], 'k','LineStyle',':','Color', [1 0 0]);
    
    % Plot wave comps
    %subplot(Nb,nsub,nsub*(n-1)+3);
    
    subplot('Position', [lmar+2*lrspc+figwid1+figwid2 yplot figwid3 fighei]);
    stem(m, abs(aD(inds)), 'Marker', '.', 'LineStyle', 'none', 'Color', [1 0 0]);

    hold on
    
    thisa = aR2*zeta(2);
    stem(m, abs(thisa(inds)), 'Marker', '.', 'Color', [0 0.5 0], 'LineStyle', '--');
    hold on
    
    thisa = aR1*zeta(1);
    stem(m, abs(thisa(inds)), 'Marker', '.', 'Color', [0 0 1]);
    
        set(gca, 'ylim', [0 0.7], 'FontSize', fontsi);
   
    
    if (n == 1)
        title('Wave Component Magnitude', 'FontSize', fontsi);
        ax = gca;
        lines = ax.Children;
        leg = legend([lines(3), lines(1), lines(2)], {'Diffracted', 'Surge', 'Heave'}, 'Position', [0.62 0.8 0.1 0.1], 'FontSize', fontsi);
        %legend('boxoff');
    end
    if (n == Nb)
        xlabel('m');
        set(gca, 'xtick', [-20 -10 0 10 20]);
    else
        set(gca, 'xtick', []);
    end

end

print('-dpng', '-r400', [folder 'pics\bristol_cyl_coefs']);

%% 4.2) Large beam range pic
load pow_big_beams

fighei = 9;

aI = exp(-1i*(-M:M)*pi/2);
aI = aI.';

h = 10;          
lam = 10;
T = IWaves.Lam2T(lam, h);
k = 2*pi/lam;

rho = 1000;
g = 9.806650;
cg = IWaves.GroupVel(T, h);
const = 2*rho*g*cg/k;

uef = IWaves.UnitEnergyFlux(rho, T, h);
uef = uef./const;

Es = zeros(NbB,1);
E2s = zeros(NbB,1);

for n = 1:NbB
    
    aS = aSsB{n}.';
    aD = 1/2*aI + aS;
    aR1 = aRsB{n,1}.';
    aR2 = aRsB{n,2}.';
    BR = [aR1 aR2];
    zeta = BR\(-aD);
    d = aD + BR*zeta;
    
    pow = (abs(aD).^2 - abs(d).^2);
    powDiff = abs(aD).^2;
    eff = pow./powDiff;
    
    E = sum(eff);
    
    Es(n) = E;
    
    totPow = sum(pow);
    CW = totPow/uef;
    
    E2 = k*CW;
    E2s(n) = E2;
end

figure;
set(gcf,'PaperPosition',[0 0 figwid fighei]);

coefs = polyfit(beamsB,Es.',1);
yfit = coefs(1)*(0:50) + coefs(2);

plot((0:50), yfit, 'k--');
hold on;
plot(beamsB, Es, 'LineWidth', 1)
axis equal

fontsi = 8;
xlabel('Beam (B/A)', 'FontSize', fontsi);
ylabel('E', 'FontSize', fontsi);
title('Power absorption (as E) v. Beam', 'FontSize', fontsi)

leg = legend('Power Absorption', 'Linear Fit');
set(leg, 'location', 'northwest', 'fontsize', fontsi)

print('-dpng', '-r400', [folder 'pics\bristol_cyl_E']);

%% 5.1) Wave field - compute

% beam = 10;
beam = 10;

M = 32;
nTheta = M*2^3;  

[aSw, aRw, ~, ~, hydBody] = wave_power_wamComp([folder 'runs'], beam, nTheta, M, 'HydroBody');

x = -30:0.25:30;
[X, Y] = meshgrid(x, x);


rho = 1000;
h = 10;          
lam = 10;
T = IWaves.Lam2T(lam, h);

M = 32;
aI = exp(-1i*(-M:M)*pi/2);

D = hydBody.DiffTM;
D = CirWaveComp.Resize2M(D{1}, M);

aD = 1/2*aI.' + aSw.';
aR1 = aRw{1};
aR2 = aRw{2};
BR = [aR1.' aR2.'];
zeta = BR\(-aD);

I = eye(2*M+1);
rhs = BR*zeta;
lhs = -(0.5*I + D);
aIR = lhs\rhs;
aSR = D*aIR;

aIR = aIR.';
aSR = aSR.';

rwaves1 = CirWaves('Out', [0 0], {aR1}, T, h);
rwaves2 = CirWaves('Out', [0 0], {aR2}, T, h);
swaves = CirWaves('Out', [0 0], {aSw}, T, h);
iwaves = CirWaves('In', [0 0], {aI}, T, h);

sRwaves = CirWaves('Out', [0 0], {aSR}, T, h);
iRwaves = CirWaves('In', [0 0], {aIR}, T, h);

rwf1 = CirWaveField(rho, rwaves1, true, X, Y);
rwf2 = CirWaveField(rho, rwaves2, true, X, Y);
swf = CirWaveField(rho, swaves, true, X, Y);
iwf = IncCirWaveField(rho, iwaves, true, X, Y);

sRwf = CirWaveField(rho, sRwaves, true, X, Y);
iRwf = IncCirWaveField(rho, iRwaves, true, X, Y);

rwf1 = rwf1*zeta(1);
rwf2 = rwf2*zeta(2);
rwf = rwf1 + rwf2;
dwf = iwf + swf;
twf = dwf + rwf;

dRwf = iRwf + sRwf;
tRwf = dRwf + rwf;

thet = linspace(0, 2*pi, 101);
cir = 5*[cos(thet).', sin(thet).'];

rwf1.RemoveGeometries(cir, 'Out');
rwf2.RemoveGeometries(cir, 'Out');
rwf.RemoveGeometries(cir, 'Out');
dwf.RemoveGeometries(cir, 'Out');
twf.RemoveGeometries(cir, 'Out');

dRwf.RemoveGeometries(cir, 'Out');
tRwf.RemoveGeometries(cir, 'Out');

etaR1 = rwf1.Elevation;
etaR2 = rwf2.Elevation;
etaR = rwf.Elevation;
etaD = dwf.Elevation;
etaT = twf.Elevation;

etaDR = dRwf.Elevation;
etaTR = tRwf.Elevation;

%% 5.2) Wave field pic

fighei = 14;

figure;
set(gcf,'PaperPosition',[0 0 figwid fighei]);

tmar = 0.05;
lmar = 0.08;

pichei = 0.215;
picwid = 0.215;

tbspc = 0.1;
lrspc = 0.05;

cols = [lmar, lmar+picwid+lrspc, lmar+2*picwid+3*lrspc];
rows = [1-tmar-pichei, 1-tmar-2*pichei-1*tbspc, 1-tmar-3*pichei-2*tbspc];

etas = [etaR1, etaR2, etaR; etaD, etaR, etaT; etaDR, etaR, etaTR];
tits = {'Surge Radiated', 'Heave Radiated', 'Total Radiated'; ...
    'Plane Diffracted', 'Total Radiated', 'Total Wave Field'; ...
    'Generalized Diffracted', 'Total Radiated', 'Total Wave Field'};

% Plot body shape    
r = 1;
beam = 10;
xs = [r, r, -r, -r, r];
ys = [-beam/2, beam/2, beam/2, -beam/2, -beam/2];

% 0 is none
% 1 if forward
% -1 is backwards
arrow1 = [-1, -1, 0; 1 0 1; 1 0 1];
arrow2 = [1 1 1; 1 1 1; 1 1 0];

for m = 1:3
    for n = 1:3
        subplot('Position', [cols(n) rows(m) picwid pichei]);
        pcolor(X, Y, real(etas{m,n}));
        
        % colormap
        thesis_cmap
        
        hold on;
        plot(cir(:,1), cir(:,2), 'k');
        set(gca, 'clim', [-1 1], 'xtick', [-20 0 20], 'ytick', [-20 0 20], 'FontSize', fontsi);
        plot(xs, ys, 'k')
        
        if ((m == 3) && (n == 1))
            xlabel('x/A', 'FontSize', fontsi);
            ylabel('y/A', 'FontSize', fontsi);
            cb = colorbar('Position', [0.915 0.2 0.03 0.6]);
            set(cb, 'FontSize', fontsi);
            title(cb, 'Re\{\eta/A\}');
        else
            set(gca,'xticklabel', [], 'yticklabel', []);
        end
        title(tits{m,n}, 'FontSize', fontsi);
        fet;
        
        dx = 0.2*picwid/2;
        len = picwid/2-3*dx;
        x1 = cols(n)+dx;
        x2 = x1+len;
        x3 = x1+picwid/2+dx;
        x4 = x3+len;
        y1 = rows(m)+pichei/2;
        y2 = y1;
        
        if (arrow1(m,n) > 0)
            annotation('arrow', 'X', [x1 x2], 'Y', [y1 y2],...
                'HeadWidth', 5, 'HeadLength', 5);
        elseif (arrow1(m, n) < 0)
            annotation('arrow', 'X', [x2 x1], 'Y', [y1 y2],...
                'HeadWidth', 5, 'HeadLength', 5);
        end
        
        if (arrow2(m,n) > 0)
            annotation('arrow', 'X', [x3 x4], 'Y', [y1 y2],...
                'HeadWidth', 5, 'HeadLength', 5);
        elseif (arrow1(m, n) < 0)
            annotation('arrow', 'X', [x4 x3], 'Y', [y1 y2],...
                'HeadWidth', 5, 'HeadLength', 5);
        end

    end
    
    
    annotation('textbox', [cols(1)+picwid+lrspc/2-0.02 rows(m)+pichei/5 0.1 0.1],...
        'String','+', 'FitBoxToText', 'on', 'LineStyle', 'none','FontSize', 20, 'FontWeight', 'bold')
    
    annotation('textbox', [cols(2)+picwid+lrspc/2+0.008 rows(m)+pichei/5 0.1 0.1],...
        'String','=','FitBoxToText', 'on', 'LineStyle', 'none','FontSize', 20, 'FontWeight', 'bold')
    
end

print('-dpng', '-r400', [folder 'pics\bristol_wfs']);

%% 6.1) Flare Radius - comp profile

lam = 10;
k = 2*pi/lam;

beam = 10;

M = 32;
nTheta = M*2^3;  

y = 0:0.1:beam/2;

volc = 1^2*pi*beam;

r = 1 + besselj(2,k*y) + besselj(3,k*y) + besselj(4,k*y)+ besselj(5,k*y);
intR2 = trapz(y,r.^2);
coef = sqrt(1/2*beam/intR2);
r = coef*r;
vol2 = 2*pi*trapz(y,r.^2);

vol2/volc

prof = [y', r'];

[aSf, aRf, compf] = wave_power_wamComp([folder 'runs'], beam, nTheta, M, 'CustProf', prof);

[aSn, aRn, compn] = wave_power_wamComp([folder 'runs'], beam, nTheta, M);

M = 32;
Mt = 20;
m = -Mt:Mt;
mn = m + 0.15*ones(size(m));
mf = m - 0.15*ones(size(m));
inds = M+1-Mt:M+1+Mt;
aI = exp(-1i*(-M:M)*pi/2);
aI = aI.';


aS = aSf.';
aDf = 1/2*aI + aS;
aR1f = aRf{1}.';
aR2f = aRf{2}.';
BRf = [aR1f aR2f];
zetaf = BRf\(-aDf);
d = aDf + BRf*zetaf;

pow = (abs(aDf).^2 - abs(d).^2);
powDiff = abs(aDf).^2;
efff = pow./powDiff;
Ef = sum(efff);

aS = aSn.';
aDn = 1/2*aI + aS;
aR1n = aRn{1}.';
aR2n = aRn{2}.';
BRn = [aR1n aR2n];
zetan = BRn\(-aDn);
d = aDn + BRn*zetan;
    
pow = (abs(aDn).^2 - abs(d).^2);
powDiff = abs(aDn).^2;
effn = pow./powDiff;
En = sum(effn);

%% 6.2) Flare radius - pic


fighei = 14;

figure;
set(gcf,'PaperPosition',[0 0 figwid fighei]);

geo = compf.Bodies.PanelGeo;
geon = compn.Bodies.PanelGeo;

subplot(2,2,1);
dep = 3;
rs = 1;
beam = 10;
h = 10;
ys = [rs, rs, -rs, -rs, +rs];
xs = [-beam/2, beam/2, beam/2, -beam/2, -beam/2];

deps = dep*ones(size(r));
plot(y, r, 'Color', [0 0 0], 'LineWidth', 1);
hold on;
plot(xs, ys, 'k--', 'LineWidth', 1);
leg = legend('Flare', 'Straight');
posL = get(leg, 'Position');
posL(1) = posL(1)-0.1;
posL(2) = posL(2)-0.25;
set(leg, 'Box', 'off', 'Position', posL, 'FontSize', 11);
plot(-beam/2*ones(size(y))+flip(y), r, 'Color', [0 0 0], 'LineWidth', 1);
plot(y, -r, 'Color', [0 0 0], 'LineWidth', 1);
plot(-beam/2*ones(size(y))+flip(y), -r, 'Color', [0 0 0], 'LineWidth', 1);
plot([-beam/2 -beam/2], [-r(end) r(end)], 'Color', [0 0 0], 'LineWidth', 1);
plot([beam/2 beam/2], [-r(end) r(end)], 'Color', [0 0 0], 'LineWidth', 1);
title('Body Profile');
axis equal
set(gca, 'ylim', [-3.5 3.5], 'xlim', [-6 6]);
axis off;


subplot(2,2,2);
% surf(geo);
% axis equal;
% set(gca, 'View', [-56 28]);
% axis off;

% Plot motions
T = IWaves.Lam2T(lam, h);
t = linspace(0,T,101);
omega = 2*pi/T;
eiomegat = exp(1i*omega*t);

xpathf = real(zetaf(1)*eiomegat);
zpathf = real(zetaf(2)*eiomegat) - dep*ones(size(eiomegat));
% arrow points
xaf = [xpathf(1), xpathf(25), xpathf(50), xpathf(75)];
zaf = [zpathf(1), zpathf(25), zpathf(50), zpathf(75)];

xpathn = real(zetan(1)*eiomegat);
zpathn = real(zetan(2)*eiomegat) - dep*ones(size(eiomegat));
% arrow points
xan = [xpathn(1), xpathn(25), xpathn(50), xpathn(75)];
zan = [zpathn(1), zpathn(25), zpathn(50), zpathn(75)];

%     delx = diff([min(xa(:)) max(xa(:))])/4;
%     dely = diff([min(za(:)) max(za(:))])/4;
%     del = delx.^2 + dely.^2;

uaf = omega*(zaf+dep*ones(size(zaf)));
vaf = -omega*xaf;
laf = sqrt(uaf.^2 + vaf.^2);
uaf = uaf./laf;
vaf = vaf./laf;

uan = omega*(zan+dep*ones(size(zan)));
van = -omega*xan;
lan = sqrt(uan.^2 + van.^2);
uan = uan./lan;
van = van./lan;
%subplot(Nb,nsub,nsub*(n-1)+2);

plot(xpathf, zpathf, 'Color', [0 0 0], 'LineWidth', 1);
hold on;
quiver(xaf, zaf, uaf, vaf, 1, 'k');

plot(xpathn, zpathn, 'k--', 'LineWidth', 1);
quiver(xan, zan, uan, van, 1, 'k');
title('Orbital Motion');

%plot(wavx, wavz, 'k--');
%plot([-30 30], [0 0], 'k:');
%plot([-30 30], [-h -h], 'k');
axis equal
set(gca, 'xlim', [-3.5 3.5], 'ylim', [-6.5 0.5]);
%set(gca, 'xlim', [-15 15], 'ylim', [-20 4]);
axis off


subplot(2,2,4);
stem(mf, efff(inds), 'Marker', '.', 'Color', [0 0 0]);
hold on;
stem(mn, effn(inds), 'Marker', '.', 'Color', [0 0.5 0], 'LineStyle', '--');
% leg = legend('Flare', 'Straight');
% posL = get(leg, 'Position');
% posL(1) = posL(1)+0.1;
% posL(2) = posL(2)+0.15;
% set(leg, 'Box', 'off', 'Position', posL);
plot([-Mt Mt], [1 1], 'k--');
% plot([E/2 E/2], [0 1], 'k','LineStyle',':','Color', [1 0 0]);
%     plot(-[E/2 E/2], [0 1], 'k','LineStyle',':','Color', [1 0 0]);
set(gca, 'ylim', [-0.1 1.1]);
title('Efficiency');
xlabel('m');
set(gca, 'xlim', [-10 10], 'xtick', [-10 -5 0 5 10]);

volf = computeVolume(geo);
voln = computeVolume(geon);

annotation('textbox', [0.36 0.54 0.1 0.1], 'String',...
    {['E_f = ' num2str(Ef, '%3.1f')], ['E_s = ' num2str(En, '%3.1f')]},...
    'LineStyle', 'none');

    
% annotation('textbox', [0.5 0.5 0.1 0.1], 'String',...
%     {['E_f = ' num2str(Ef, '%3.1f'), '      E_n = ' num2str(En, '%3.1f')],...
%     ['Vol_f = ' num2str(volf, '%3.1f'),'      Vol_n= ' num2str(voln, '%3.1f')],...
%     ['E_f/Vol_f = ' num2str(Ef/volf, '%5.3f'), '      E_n/Vol_n= ' num2str(En/voln, '%5.3f')]},...
%     'LineStyle', 'none');

 subplot(2,2,3);
 %stem(m, abs(aDf(inds)), 'Marker', '.', 'LineStyle', 'none', 'Color', [1 0 0]);
 hold on
 thisa = aR2f*zetaf(2);
 stem(mf, abs(thisa(inds)), 'Marker', '.', 'Color', [0 0 0]);
 hold on
 thisa = aR1f*zetaf(1);
 stem(mf, abs(thisa(inds)), 'Marker', '.', 'Color', [0 0 0]);
 
 thisa = aR1n*zetan(2);
 stem(mn, abs(thisa(inds)), 'Marker', '.', 'Color', [0 0.5 0], 'LineStyle', '--');
 hold on
 thisa = aR2n*zetan(1);
 stem(mn, abs(thisa(inds)), 'Marker', '.', 'Color', [0 0.5 0], 'LineStyle', '--');
 set(gca, 'xlim', [-10 10], 'ylim', [0 1], 'Box', 'on');
 
title('Wave Component Magnitude');
xlabel('m');
set(gca, 'xlim', [-10 10], 'xtick', [-10 -5 0 5 10]);
 
 plot([-Mt Mt], [0.5 0.5], 'k--');
 
print('-dpng', '-r400', [folder 'pics\bristol_flare']);