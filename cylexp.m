%% Cylindrical Experiment
% copied from ewtec_PaperPics in
% C:\Users\s1213969\Dropbox\matlab\studies\exp_res\expRes_ewtec_paper

folder = 'C:\Users\s1213969\Dropbox\matlab\thesis\';
dataDir = 'N:\RDS\Cylindrical_wave_field_experiments\Raw_data\';

fontsi = 10;
figwid = 14;

bodies = {'Flap', 'Atten'};
waveTypes = {'Scat', 'Rad'};
orients = [0 45];
fs = [0.8 1 1.25];
M = 5;

%% Evanescent wave check

load cylexp_cyl_grid;
load cylexp_cyl_hb;

T = cylComp.T;
h = cylComp.H;
cylWF.BodyMotions = cylComp.Motions;

[X, Y] = cylWF.FieldPoints;

iwaves(1) = PlaneWaves(1, T, 0, h);
iwaves(2) = PlaneWaves(1, T, 7/4*pi, h);

cylComp2 = HydroBodyComp(cylHydBody, iwaves);
cylWF2 = cylComp2.WaveField(true, X, Y, 'NoVel');

etaWc = cylWF.Elevation('Diffracted');
etaCc = cylWF2.Elevation('Diffracted');

load cylexp_flap_grid;
load cylexp_flap_hb;

flapWF.BodyMotions = flapComp.Motions;

flapComp2 = HydroBodyComp(flapHydBody, iwaves);
flapWF2 = flapComp2.WaveField(true, X, Y, 'NoVel');

etaWf = flapWF.Elevation('Diffracted');
etaCf = flapWF2.Elevation('Diffracted');

delC = cell(3,1);
delF = cell(3,1);
for m = 1:3
    delC{m} = 100*abs(etaCc{m,1}-etaWc{m,1});
    delF{m} = 100*abs(etaCf{m,1}-etaWf{m,1});
end

%%

figure;
set(gcf,'PaperPosition',[0 0 12 14]);

lmar = 0.15;
tmar = 0.06;
picwid = 0.3;
pichei = 0.25;
lrspc = 0.05;
tbspc = 0.05;

r = 0.8;
ctheta = linspace(0,2*pi,101);

thet = linspace(-pi/2,pi/2,31);
curv = [0.4*ones(size(thet)); zeros(size(thet))] + 0.08*[cos(thet); sin(thet)];
bodyC = [-curv, curv, -curv(:,1)];
bodyC = bodyC';

bodyF = [-0.04 0.3; 0.04 0.3; 0.04 -0.3; -0.04 -0.3; -0.04 0.3];

for m = 1:3
    ypos = 1-tmar-m*pichei-(m-1)*tbspc;
    %subplot(3,2,(m-1)*2+1);
    subplot('Position', [lmar ypos picwid pichei]);
    pcolor(X,Y,delC{m});
    fet;
    hold on;
    plot(bodyC(:,1), bodyC(:,2), 'k', 'LineWidth', 1); axis equal
    plot(r*cos(ctheta), r*sin(ctheta), 'w:', 'LineWidth', 1);
    set(gca, 'clim', [0 5], 'FontSize', fontsi, 'xtick', [], 'ytick', []);
    if (m == 1)
        title('Attenuator', 'FontSize', fontsi);
    end
    if (m == 3)
        set(gca, 'xtick', [-2 -1 0 1 2], 'ytick', [-2 -1 0 1 2])
        ylabel({[num2str(1./T(m)) ' Hz'], ['y (m)']}, 'FontSize', fontsi);
        xlabel('x (m)', 'FontSize', fontsi);
    else
        ylabel([num2str(1./T(m)) ' Hz'], 'FontSize', fontsi);
    end

    %subplot(3,2,(m-1)*2+2);
    subplot('Position', [lmar+picwid+lrspc ypos picwid pichei]);
    pcolor(X,Y,delF{m});
    fet;
    hold on;
    plot(bodyF(:,1), bodyF(:,2), 'k', 'LineWidth', 1); axis equal
    plot(r*cos(ctheta), r*sin(ctheta), 'w:', 'LineWidth', 1);
    set(gca, 'clim', [0 5], 'FontSize', fontsi, 'xtick', [], 'ytick', []);
    if (m == 1)
        title('Flap', 'FontSize', fontsi);
    end
end

% WF colorbar

caxis = colorbar('location','manual','position',[0.87 0.3 0.04 0.4]);
set(caxis, 'fontsize', fontsi);
%xlabel(caxis, '|\eta/a|', 'fontsize', fontsi);
x1 = caxis.Position(1);
y1 = caxis.Position(2);
dxcl = -0.06;
dycl = -0.08;
annotation(gcf, 'textbox', [x1+dxcl+0.02 y1+dycl+0.01 0.2 0.07], 'string', {'Err (%)'}, 'linestyle', 'none', 'fontsize', fontsi);

thesis_cmap

print('-dpng','-r500', [folder 'pics\evan_err'])

%% Compute result
%{
body = 'Flap';
body = 'Atten';
body = 'None'

motion = 'Fix';
motion = 'Rad';
motion = 'Free';

orient = 0;
orient = 22.5;
orient = 45;
orient = 67.5;
orient = 90;
orient = -22.5;
orient = -45;
orient = -67.5;

Optional arguments
'Wave'
'Motion'
'WGSetup'
'Orient'
'StartT'
%}

plotChans = {};

ignoreRuns = [];

body = bodies{2};
waveType = waveTypes{2};
fR = fs(3);
orient = orients(1);

if (strcmp(waveType, 'Scat'))
    [Ameas, Acyl, Agrid, Awam, am, amW, X, Y, R2c, R2w, Std, incAmp, motAmp] = ...
        cylexp_analyzeData(dataDir, body, fR, waveType, 'Orient', orient, 'M', M);
    amp = incAmp;
else
    [Ameas, Acyl, Agrid, Awam, am, amW, X, Y, R2c, R2w, Std, incAmp, motAmp] = ...
        cylexp_analyzeData(dataDir, body, fR, waveType, 'M', M);
    amp = motAmp;
end

% Make a figure

cylexp_makePlot(body, waveType, fR, Ameas, Acyl, am{1}, Agrid, X, Y,...
    'PlotC', 'Orient', orient, 'Modelam', amW, 'ModelA', Awam, 'R2', R2c, R2w, 'IncMot', amp)

fstr = num2str(fR);
idot = strfind(fstr,'.');
if (~isempty(idot))
    fstr = fstr([1:idot-1,idot+1:end]);
end

fileName = [folder 'pics\' body '_' waveType '_' fstr '_' num2str(orient) '_m'];

print('-dpng','-r400', fileName)


%% phase picture

body = 'Flap';
fR = 1;
waveType = 'Rad';

[Ameas, Acyl, ~, Awam] = cylexp_analyzeData(dataDir, body, fR, waveType, 'M', M);

sind = 1;
r0 = 60:20:200;

figure;
set(gcf,'PaperPosition',[0 0 15 10]);

lmar = 0.12;
tmar = 0.1;
figwid = 0.85;
fighei = 0.32;
tbspc = 0.1;

subplot('Position', [lmar 1-tmar-fighei figwid fighei]);

plot(r0, real(Ameas{sind}), 'LineStyle','--', 'Marker', 'x', 'Color', [0 0 0], 'DisplayName', 'Measured');
hold on;
plot(r0, real(Acyl{sind}), 'LineStyle','-.', 'Marker', '+', 'Color', [0 0 1], 'DisplayName', 'Cylindrical');
plot(r0, real(Awam{sind}), 'LineStyle',':', 'Marker', '*', 'Color', [1 0 0], 'DisplayName', 'WAMIT');

set(gca, 'ylim', [-3 3], 'xlim', [50 210], 'xticklabel', [], 'FontSize', fontsi);
title('Flap, Radiated, 1.00 Hz - s1','FontSize', fontsi);
ylabel('Re\{\eta\} (cm)','FontSize', fontsi);

subplot('Position', [lmar 1-tmar-2*fighei-tbspc figwid fighei]);

plot(r0, imag(Ameas{sind}), 'LineStyle','--', 'Marker', 'x', 'Color', [0 0 0], 'DisplayName', 'Measured');
hold on;
plot(r0, imag(Acyl{sind}), 'LineStyle','-.', 'Marker', '+', 'Color', [0 0 1], 'DisplayName', 'Cylindrical');
plot(r0, imag(Awam{sind}), 'LineStyle',':', 'Marker', '*', 'Color', [1 0 0], 'DisplayName', 'WAMIT');
set(gca, 'ylim', [-3 3], 'xlim', [50 210],'FontSize', fontsi);
ylabel('Im\{\eta\} (cm)','FontSize', fontsi);
xlabel('radial position (cm)','FontSize', fontsi);
leg = legend('Show');
set(leg, 'Position', [0.73 0.48 0.2 0.15],'FontSize', fontsi);

print('-dpng','-r300', [folder 'pics\phase_s' num2str(sind) '_' body '_' waveType '_f100_' num2str(orient) '_m'])

%% Get data for Fourier trans pics

fstr = 'f100';

runs = ['Test_Fl_Fix_b0_' fstr '_pw_wg2_*.csv'];
names = dir([dataDir runs]);
Fnames = {names.name};

runs = ['Test_Fl_Rad_' fstr '_wg2_*.csv'];
names = dir([dataDir runs]);
Rnamez = {names.name};

trun = 2;

nR = 1;
nF = 2;
o = 13;

data = wfe_load_data_file([dataDir Rnamez{nR}]);   % load the data file

cDataR = wfe_get_cald_data(dataDir, data);      % calibrate data
sampleFreq = data.SampleFreq;       
startT = 40;

tR = data.time;

[fR, aR, epR, offR, SR, freqsR, a2R, ep2R] = wfe_compute_sig_vals(sampleFreq, cDataR(:,o), 'StartTime', startT, 'Range');

data = wfe_load_data_file([dataDir Fnames{nF}]);   % load the data file
data.WGNames{o}

cDataF = wfe_get_cald_data(dataDir, data);      % calibrate data
sampleFreq = data.SampleFreq;       

[fF, aF, epF, offF, SF, freqsF, a2F, ep2F] = wfe_compute_sig_vals(sampleFreq, cDataF(:,o), 'StartTime', startT);

tF = data.time;

%% Plot Fourier trans pics

figure;
set(gcf,'PaperPosition',[0 0 12 14]);

tmar = 0.06;
lmar = 0.16;

fighei = 0.18;
figwid = 0.34;

tbspc = 0.1;
lrspc = 0.1;

for n = 1:2
    if (n == 1)
        origSig = cDataF(:,o).';
        t = tF;
        S = SF;
        freqs = freqsF;
        f = fF;
        a = aF;
        ep = epF;
        a2 = a2F;
        ep2 = ep2F;
        row = 0;
        
        figPos = [lmar 1-tmar-fighei figwid fighei];
    else
        origSig = cDataR(:,o).';
        t = tR;
        S = SR;
        freqs = freqsR;
        f = fR;
        a = aR;
        ep = epR;
        a2 = a2R;
        ep2 = ep2R;
        row = 3;
        
        figPos = [lmar+lrspc+figwid 1-tmar-fighei figwid fighei];
    end
    
    subplot('Position', figPos);
    plot(t, origSig, 'LineWidth', 0.5);
    hold on;
    [iUp, iDown] = meanCrossPeaks(origSig);

%     plot(t(iUp), origSig(iUp), 'k', 'LineWidth', 0.5);
%     plot(t(iDown), origSig(iDown), 'k', 'LineWidth', 0.5);
    set(gca, 'FontSize', fontsi, 'xlim', [0 150]);
    if (n == 2)
        %set(gca, 'yticklabel', []);
        title('Radiated (c12)', 'FontSize', fontsi)
    else
        ylabel({'Total Time', 'Elevation (cm)'}, 'FontSize', fontsi);
        title('Diffracted (c12)', 'FontSize', fontsi)
    end
    xlabel('Time (s)', 'FontSize', fontsi);
    
    % plot Spectrum
    figPos(2) = figPos(2) - tbspc - fighei;
    subplot('Position', figPos);
    plot(freqs, abs(S), 'LineWidth', 0.5);

    set(gca, 'FontSize', fontsi, 'xlim', [0 5], 'xtick', (0:5));
    if (n == 2)
        %set(gca, 'yticklabel', []);
    else
        ylabel({'Spectrum', 'Amplitude (cm)'}, 'FontSize', fontsi);
    end
    xlabel('Frequency (Hz)', 'FontSize', fontsi);

    % plot close up of signal with linear and second order fit
    linSig = a*cos(2*pi*f*t + ep);
    secSig = linSig + a2*cos(2*2*pi*f*t + ep2);

    % just 2 seconds in the middle of the run
    %startT2 = (t(end) - startT)/2;
    startT2 = startT;
    [~, startI] = min(abs(t - startT2));
    stopT = startT2 + trun; 
    if (stopT > t(end))
        stopT = t(end);
    end
    [~, stopI] = min(abs(t - stopT));

     figPos(2) = figPos(2) - tbspc - fighei;
    subplot('Position', figPos);
    plot(t(startI:stopI), origSig(startI:stopI), 'Linewidth', 1);
    hold on;
    plot(t(startI:stopI), linSig(startI:stopI), 'Color', [0 0.5 0],'LineStyle','-.', 'LineWidth', 0.5);
    plot(t(startI:stopI), secSig(startI:stopI), 'k--', 'LineWidth', 0.5);

    leg1 = legend({'Data', 'Linear', '2nd Order'});
    set(leg1,'Position',[0.52 0.03 0.07 0.072], 'FontSize', fontsi);
    set(gca, 'FontSize', fontsi, 'xlim', [t(startI) t(stopI)]);
    if (n == 2)
        %set(gca, 'yticklabel', []);
    else
        ylabel({'Time', 'Elevation (cm)'}, 'FontSize', fontsi);
    end
    xlabel('Time (s)', 'FontSize', fontsi);
end

print('-dpng','-r300', [folder 'pics\time_fourier_pic'])

























