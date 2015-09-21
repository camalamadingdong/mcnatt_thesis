function [] = cylexp_makePlot(body, wave, f, Ameas, Acyl, am, Agrid, X, Y, varargin)

[opts, args] = checkOptions({{'PlotC'}, {'Orient', 1}, {'Modelam', 1}, ...
    {'ModelA', 1}, {'ShowReal'}, {'R2', 2}, {'Std', 1}, {'IncMot',1}}, varargin);

if (opts(1))
    nrow = 6;
    rowSt = 1;
else
    nrow = 5;
    rowSt = 0;
end

if (opts(2))
    orient = args{2};
else
    orient = 0;
end

if (opts(3))
    amW = args{3};
else
    amW = 0;
end

if (opts(4))
    Awam = args{4};
else
    Awam = [];
end

if (opts(6))
    R2c = args{6}{1};
    R2w = args{6}{2};
else
    R2c = 0;
    R2w = 0;
end

if (opts(7))
    std = args{7};
else
    std = [];
end

if (opts(8))
    motAmp = args{8};
else
    motAmp = [];
end

figure;
set(gcf,'PaperPosition',[0 0 14 17]);
fntsize = 10;
%set(gcf,'PaperPosition',[0 0 18 12]);

tmar = 0.05;
lmar = 0.12;

fighei1 = 0.12;
fighei2 = 0.14;
wfhei = 0.45;
wfwid = 0.5;
figwid0 = 0.3;
figwid1 = 0.7;
figwid2 = 0.2;

tbspc1 = 0.015;
tbspc2 = 0.09;
tbspc3 = 0.13;
tbspc4 = 0.11;
lrspc = 0.08;


maxA = 0;
for n = 1:nrow
    mAn = max(abs(Ameas{n}));
    if (mAn > maxA)
        maxA = mAn;
    end
    
    mAn = max(abs(Acyl{n}));
    if (mAn > maxA)
        maxA = mAn;
    end
    
    mAw = max(abs(Awam{n}));
    if (mAw > maxA)
        maxA = mAw;
    end
end

spokes = {'s1', 's2', 's3', 's4', 's5'};
units = 100; % 100 - cm, 1 - m, 1000 - mm
r = units*(0.6:0.2:2);
xlims = [r(1) r(end)] + units*[-0.1 0.1];

if (strcmp(body, 'Atten'))
    body = 'Attenuator';
end

if (strcmp(wave, 'Rad'))
    tit = [body ', Radiated, ' num2str(f,'%4.2f') ' Hz'];
    if (~isempty(motAmp))
        tit = [tit ', |\xi| = ' num2str(motAmp,'%3.1f') '^{\circ}'];
    end
else
    tit = [body ', ' num2str(orient) ' deg, Scattered, ' num2str(f,'%4.2f') ' Hz'];
    if (~isempty(motAmp))
        tit = [tit ', |\eta^{I}| = ' num2str(motAmp,'%3.1f') ' cm'];
    end
end

if (opts(5))
    plotfun = @real;
    ylim = maxA*1.1*[-1 1];
else
    plotfun = @abs;
    ylim = maxA*1.1*[0 1];
end

sRad = [0, pi, 5*pi/4, 3*pi/2, 7*pi/4]; 

if (rowSt == 1)
    %subplot(nrow,2,1:2);
    subplot('Position', [lmar 1-tmar-fighei1 figwid1 fighei1]);
    thet = linspace(0,2*pi*(1-1/24),24);
    if (~isempty(std))
        errorbar(thet, plotfun(Ameas{6}), abs(std{6}), 'LineStyle','--', 'Marker', 'x', 'Color', [0 0 0], 'DisplayName', 'Measured');
    else
        plot(thet, plotfun(Ameas{6}), 'LineStyle','--', 'Marker', 'x', 'Color', [0 0 0], 'DisplayName', 'Measured');
    end
    hold on
    plot(thet, plotfun(Acyl{6}), 'LineStyle','-.', 'Marker', '+', 'Color', [0 0 1], 'DisplayName', 'Cylindrical');
    if (~isempty(Awam))
        plot(thet, plotfun(Awam{6}), 'LineStyle',':', 'Marker', '*', 'Color', [1 0 0], 'DisplayName', 'WAMIT');
    end
        
    set(gca, 'xlim', [-pi/12 2*pi], 'ylim', ylim, 'FontSize', fntsize);
    title(tit, 'FontSize', fntsize);
    xlabel('Angular position (rad)', 'FontSize', fntsize);
    ylabel({'c', 'amp (cm)'}, 'FontSize', fntsize);
    leg1 = legend(gca,'show');
    set(leg1,'Position',[0.83 0.86 0.114 0.072], 'FontSize', fntsize);
    for m = 1:5
        plot([sRad(m) sRad(m)], ylim, 'LineStyle','--', 'LineWidth', 0.2, 'Color', [0 0 0]);
    end
end

figpos = [lmar 1-tmar-2*fighei1-tbspc2 figwid0 fighei1];

for n = 1:5
    %subplot(nrow,2,2*(n-1+rowSt)+1);
    
    subplot('Position', figpos);
    figpos(2) = figpos(2) - tbspc1 - fighei1;
    
    if (~isempty(std))
        errorbar(r, plotfun(Ameas{n}), abs(std{n}), 'LineStyle','--', 'Marker', 'x', 'Color', [0 0 0]);
    else
        plot(r, plotfun(Ameas{n}), 'LineStyle','--', 'Marker', 'x', 'Color', [0 0 0]);
    end
    hold on
    plot(r, plotfun(Acyl{n}), 'LineStyle','-.', 'Marker', '+', 'Color', [0 0 1]);
    if (~isempty(Awam))
        plot(r, plotfun(Awam{n}), 'LineStyle',':', 'Marker', '*', 'Color', [1 0 0]);
    end
    
    plot([r(2) r(2)], ylim, 'LineStyle','--', 'LineWidth', 0.2, 'Color', [0 0 0]);
    set(gca, 'xlim', xlims, 'ylim', ylim, 'FontSize', fntsize);

    ylabel({spokes{n}, 'amp (cm)'}, 'FontSize', fntsize);
    if (n == 1)
        if (rowSt == 0)
            title(tit, 'FontSize', fntsize)
        end
        if (R2c ~= 0)
            annotation('textbox', [0.44,0.1,0.1,0.1],...
           'String', {['R^2_c = ' num2str(R2c, '%4.2f')], ['R^2_w = ' num2str(R2w, '%4.2f')]}, 'FontSize', fntsize);
        end
    end
    
    if (n == 5)
        xlabel('radial position (cm)', 'FontSize', fntsize);
    else
        set(gca, 'xticklabel', [], 'FontSize', fntsize);
    end
end
        
M = length(am);
M = (M - 1)/2;

%subplot(nrow, 2, 2*(rowSt+1));
figpos = [lmar+figwid0+lrspc 1-tmar-2*fighei1-tbspc3 figwid2 fighei2];
subplot('Position', figpos);
stem((-M:M) - 0.1*ones(1,2*M+1), abs(am), 'Marker', '.','LineStyle','-.','Color', [0 0 1],'LineWidth', 1)
if (length(amW) > 1)
    hold on;
    stem((-M:M) + 0.1*ones(1,2*M+1), abs(amW), 'Marker', '.','LineStyle',':','Color', [1 0 0],'LineWidth', 1)
end
set(gca, 'xlim', [-M M], 'xtick', (-5:5), 'xticklabel', {[], '-4', [], '-2', [], '0', [], '2', [], '4', []}, 'FontSize', fntsize)
xlabel('m', 'FontSize', fntsize);
ylabel('(cm)', 'FontSize', fntsize);
title('|b_m|', 'FontSize', fntsize);

theta = linspace(-pi, pi, 101);
afun = zeros(size(theta));
afunW = zeros(size(theta));

for n = -M:M
    afun = afun + exp(1i*n*pi/2)*am(n+M+1)*exp(1i*n*theta);
    if (length(amW) > 1)
        afunW = afunW + exp(1i*n*pi/2)*amW(n+M+1)*exp(1i*n*theta);
    end
end

afun = 1/pi*afun;
afunW = 1/pi*afunW;

%subplot(nrow, 2, 2*(rowSt+2));
figpos(1) = figpos(1) + lrspc + figwid2;
subplot('Position', figpos);
plot(theta, abs(afun),'LineStyle','-.','Color', [0 0 1], 'LineWidth', 1);
if (length(amW) > 1)
    hold on;
    plot(theta, abs(afunW),'LineStyle',':','Color', [1 0 0],'LineWidth', 1);
end
set(gca, 'xlim', [-pi pi], 'FontSize', fntsize)
ylabel('(cm)', 'FontSize', fntsize);
xlabel('Direction (rad)', 'FontSize', fntsize)
title('|F(\theta)|', 'FontSize', fntsize);

%subplot(nrow,2,2*(rowSt*[1,1,1] + [3:5]));
figpos(1) = figpos(1) - lrspc - figwid2 + 0.03;
figpos(2) = figpos(2) - tbspc4 - wfhei;
figpos(3:4) = [wfwid wfhei];
subplot('Position', figpos);

pcolor(X, Y, real(Agrid));
thesis_cmap
cbax = colorbar;
xlabel(cbax, 'Re\{\eta\} (cm)', 'FontSize', fntsize);
axis off;
set(gca, 'clim', ylim(2)*[-1 1], 'FontSize', fntsize);
title('Wave field from measured coefficients (t=0)', 'FontSize', fntsize);
fet;
hold on;
tank = wfe_tank_points();
plot(tank(:,1), tank(:,2), 'k', 'LineWidth', 1); axis equal
if (strcmp(body, 'Attenuator'))
    thet = linspace(-pi/2,pi/2,31);
    curv = [0.4*ones(size(thet)); zeros(size(thet))] + 0.08*[cos(thet); sin(thet)];
    bodyPnts = [-curv, curv, -curv(:,1)];
    bodyPnts = bodyPnts';
else
    bodyPnts = [-0.04 0.3; 0.04 0.3; 0.04 -0.3; -0.04 -0.3; -0.04 0.3];
end
if (orient ~= 0)
    R = [cosd(orient) -sind(orient); sind(orient) cosd(orient)];
    N = size(bodyPnts, 1);
    for n = 1:N
        newPnts = R*bodyPnts(n,:)';
        bodyPnts(n,:) = newPnts';
    end
end
plot(bodyPnts(:,1), bodyPnts(:,2), 'k', 'LineWidth', 1); axis equal


cylexp_plotWg(2, 'w', 'LineWidth', 1);

x1 = 0.59;
y1 = 0.22;
w1 = 0.1;
h1 = 0.1;
annotation('textbox', [x1, y1, w1, h1], 'String', 's1', 'LineStyle', 'none', 'FontSize', fntsize);
annotation('textbox', [x1+0.01, y1-0.042, w1, h1], 'String', 's2', 'LineStyle', 'none', 'FontSize', fntsize);
annotation('textbox', [x1+0.075, y1-0.066, w1, h1], 'String', 's3', 'LineStyle', 'none', 'FontSize', fntsize);
annotation('textbox', [x1+0.13, y1-0.042, w1, h1], 'String', 's4', 'LineStyle', 'none', 'FontSize', fntsize);
annotation('textbox', [x1+0.17, y1, w1, h1], 'String', 's5', 'LineStyle', 'none', 'FontSize', fntsize);

end