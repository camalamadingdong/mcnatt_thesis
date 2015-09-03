%% Introduction pictures
% copied from jfm_PaperPics in N:\RDS\WAMIT\array_study\jfm
% 0) set up stuff

folder = 'C:\Users\s1213969\Dropbox\matlab\thesis\';
fontsi = 10;
figwid = 7;
fighei = 5;

% sine wave
a = 1;
lam = 10;
k = 2*pi/lam;
h = 20;
T = IWaves.Lam2T(lam, h);
beta = 0;
rho = 1000;

x = linspace(0, 30, 301);

wav2d = a*cos(k*x);

figure;
set(gcf, 'PaperPosition', [0 0 figwid fighei]);
plot(x,wav2d,'Linewidth',2)
fet;
axis off
print('-dpng', '-r400', [folder 'pics\wav2d']);

% plane wave
y = linspace(-20, 20, 401);
[X, Y] = meshgrid(x, y);

wavp = PlaneWaves(a, T, beta, h);
pwf = PlaneWaveField(rho, wavp, true, X, Y);

eta = pwf.Elevation;

figure;
set(gcf, 'PaperPosition', [0 0 figwid fighei]);
surf(X,Y,real(eta{1}));
fet;
axis off
print('-dpng', '-r400', [folder 'pics\plane_wave']);

% cylindrical wave
M = 6;
s = 20;
k = 2*pi/lam;

am = cirWaveCoefsFromSpread(M, s, k);
am = 6*am;

wavc = CirWaves('In', [0 0], {am}, T, h);

x = linspace(0, 60, 301);
y = linspace(-50, 50, 501);
[X, Y] = meshgrid(x, y);

cwf = IncCirWaveField(rho, wavc, true, X, Y);

eta = cwf.Elevation;

figure;
set(gcf, 'PaperPosition', [0 0 figwid fighei]);
surf(X,Y,real(eta{1}));
fet;
set(gca, 'Color', 'none');
axis off
print('-dpng', '-r400', [folder 'pics\curve_wave']);
