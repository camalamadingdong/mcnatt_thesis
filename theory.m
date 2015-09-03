%% Theory 

folder = 'C:\Users\s1213969\Dropbox\matlab\thesis\pics\';

fontsi = 10;
figwid = 14;

%% graph of Bessel functions

m1 = [0 1];
m2 = 7;
kr1 = linspace(0,15,101);
kr2 = linspace(0,6,101);

[M1, KR] = meshgrid(m1, kr1);

J1 = besselj(M1,KR);
Y1 = bessely(M1,KR);
J2 = besselj(m2, kr1);
Y2 = bessely(m2, kr1);

[M1, KR2] = meshgrid(m1, kr2);

K1 = besselk(M1, KR2);
I1 = besseli(M1, KR2);
K2 = besselk(m2, kr2);
I2 = besseli(m2, kr2);

figure;
set(gcf,'PaperPosition',[0 0 15 15]);

subplot(2,2,1);
l0 = plot(kr1,J1);
hold on;
l2 = plot(kr1,Y1,'--');
set(l2(1), 'Color', l0(1).Color)
set(l2(2), 'Color', l0(2).Color)
set(gca, 'ylim', 2*[-1 1], 'ytick', [-2 -1 0 1 2]);
leg = legend('J_0', 'J_1', 'Y_0', 'Y_1');
set(leg, 'Box', 'off', 'Position', leg.Position + [0 0.054 0 0]);
%ylabel('Magnitude');
xlabel('kr');

subplot(2,2,2);
l1 = plot(kr1,J2);
hold on;
l2 = plot(kr1,Y2,'--');
set(l2(1), 'Color', l0(1).Color)
set(gca, 'ylim', 2*[-1 1], 'ytick', [-2 -1 0 1 2]);
leg = legend('J_7', 'Y_7');
set(leg, 'Box', 'off');
%ylabel('Magnitude');
xlabel('kr');

subplot(2,2,3);
l0 = plot(kr2,K1);
hold on;
l2 = plot(kr2,I1,'--');
set(l2(1), 'Color', l0(1).Color)
set(l2(2), 'Color', l0(2).Color)
set(gca, 'ylim', 2*[-0.1 1], 'ytick', [0 1 2]);
leg = legend('K_0', 'K_1', 'I_0', 'I_1');
set(leg, 'Box', 'off');
%ylabel('Magnitude');
xlabel('kr');

subplot(2,2,4);
l1 = plot(kr2,K2);
hold on;
l2 = plot(kr2,I2,'--');
set(l2(1), 'Color', l0(1).Color)
set(gca, 'ylim', 2*[-0.1 1], 'ytick', [0 1 2]);
leg = legend('K_7', 'I_7');
set(leg, 'Box', 'off', 'Location', 'Northwest');
%ylabel('Magnitude');
xlabel('kr');

print('-dpng','-r300', 'C:\Users\s1213969\Dropbox\matlab\thesis\pics\bessels')

%% Plot of outgoing partial waves

fighei = 15;

r0 = 2*sqrt(2);
vang = [-35.5 54];
xl = [-r0 r0];
yl = xl;
zl = [-1.5 1.5];
cl = [-0.8 0.8];

loc = [0 0];
rho = 1;

x = -2:0.01:2;
y = x;
[X Y] = meshgrid(x, y);

theta = 0:2*pi/100:2*pi;
xt = r0*cos(theta);
yt = r0*sin(theta);

lam = 1;
h = 10;
T = IWaves.Lam2T(lam, h);

Bs = {1, [-0.5 0 0.5], [0.5 0 0 0 0.5], [-0.5 0 0 0 0 0 0.5]};

figure;

set(gcf, 'PaperPosition', [0 0 figwid fighei]);

tmar = 0.1;
lmar = 0.1;

wid = (1-lmar)/2;
hei = (1-tmar)/4;

lefs = lmar:wid:1;
bots = (1-tmar-hei):-hei:0;

xa = tmar/2-0.02;
dya = hei/2-0.05;

for m = 1:4
    


    B = Bs{m};
    wc = CirWaves('Out', [0 0], {B}, T, h);

    w = CirWaveField(rho, wc, true, X, Y);
    eta = w.Elevation;


    for n = 1:2
        subplot('Position', [lefs(n), bots(m), wid, hei]);
        if (n == 1)
            surf(X, Y, real(eta{1}));
        else
            surf(X, Y, imag(eta{1}));
        end
        shading flat;
        view(gca, vang);
        set(gca, 'zlim', zl, 'xlim', xl, 'ylim', yl, 'xtick', [], 'ytick', [], 'ztick', [])
        set(gca, 'clim', cl)

        f = cos((m-1)*theta);

        hold on;
        plot3(xt, yt, f);

        axis off
    end
    
    if (m == 1)
        mstr = ['m = ', num2str(m-1)];
    else
        mstr = ['m = |', num2str(m-1) '|'];
    end
    annotation(gcf, 'textbox', [xa bots(m)+dya 0.1 0.1], 'string', mstr, 'linestyle', 'none', 'fontsize', fontsi)
    
end

annotation(gcf, 'textbox', [lmar+wid/2-0.05 1-tmar-0.05 0.1 0.1], 'string', 'Re\{\eta\}', 'linestyle', 'none', 'fontweight', 'bold', 'fontsize', fontsi)
annotation(gcf, 'textbox', [lmar+wid+wid/2-0.05 1-tmar-0.05 0.1 0.1], 'string', 'Im\{\eta\}', 'linestyle', 'none', 'fontweight', 'bold', 'fontsize', fontsi)


%print('-dpng', '-r300', [path '\partWaves']);
print('-dpng', '-r400', [folder 'partOutWaves']);

%% Plot of incident partial waves

fighei = 5;

figure; 
set(gcf, 'PaperPosition', [0 0 figwid fighei]);

tmar = 0.25;
phei = 0.7;
pwid = 0.17;
lmar = 0.025;
lrspc = 0.025;

lam = 1;
h = 1;
k = 2*pi/lam;

M1 = 8;

beta = 0;
a0 = IHydroComp.IAmps(M1, [0 0], k, beta);

x = -2*lam:lam/100:2*lam;
y = x;

[X, Y] = meshgrid(x, y);
Z = zeros(size(X));

waveT = zeros(size(X));

r = sqrt(X.^2 + Y.^2);
theta = atan2(Y, X);

for m = 0:M1
    if (m == 0)
        psi = besselj(m, k*r).*exp(1i*m*theta);

        wavem = a0(M1+1+m)*psi;
    else
        psip = besselj(m, k*r).*exp(1i*m*theta);
        psin = besselj(-m, k*r).*exp(-1i*m*theta);
        
        wavemp = a0(M1+1+m)*psip;
        wavemn = a0(M1+1-m)*psin;
        wavem = wavemp + wavemn;
    end
    
    waveT = waveT + wavem;
    
    if (mod(m,2) == 0)
        subplot('Position', [lmar+(m/2)*(pwid+lrspc) 1-tmar-phei pwid phei]);
        pcolor(X, Y, real(waveT));
        fet;
        shading interp;
        set(gca, 'xtick', [], 'ytick', [], 'clim', [-1 1]);
        if (m == 4)
            title({'Re\{\eta\}', '',['M = ' num2str(m)]},'FontSize', fontsi);
        else
            title(['M = ' num2str(m)],'FontSize', fontsi);
        end
    end
end

print('-dpng', '-r400', [folder '\incCylWaves']);