function [Mrec] = cylcomp_checkEvCut(waveC, waveG, wec, varargin )

if (~isempty(varargin))
    clim = varargin{1};
else
    clim = [0 0.02];
end

M = 15;
N = 20; 

[X Y] = waveG.FieldPoints;
R = sqrt(X.^2 + Y.^2);
Err = zeros(1,M+1);

waveG.BodyMotions = wec.Motions;
eta = waveG.Elevation('Total');
etaG = eta{1};

figure;

% M = 0
m = 0;
[waveC2 waveF r0] = cylcomp_makeCirWF(waveC, m, N, X, Y);
thetar = 0:2*pi/100:2*pi;
cirx = r0*cos(thetar);
ciry = r0*sin(thetar);

inds = (R > r0);
Npts = sum(sum(inds));
waveC2.BodyMotions = wec.Motions;

eta = waveC2.Elevation('Total');
etaC = eta{1};

err = abs(etaC - etaG)./abs(etaG);
inan = find(isnan(err));
err(inan) = 0;
Err(1) = sum(err(inds));

subplot(3,3,1);
pcolor(X,Y,err);
shading interp;
axis equal;
axis tight;
set(gca, 'clim', clim)

hold on;
plot(cirx, ciry, 'w');

for m = 0:M

    waveC2 = cylcomp_makeCirWF(waveC, m, N, X, Y);
    waveC2.BodyMotions = wec.Motions;
    
    eta = waveC2.Elevation('Total');
    etaC = eta{1};
    
    err = abs(etaC - etaG)./abs(etaG)./Npts;
    inan = find(isnan(err));
    err(inan) = 0;
    Err(1+m) = sum(err(inds));

    subplot(4,4,m+1);
    pcolor(X,Y,100*err);
    shading interp;
    axis equal;
    axis tight;
    set(gca, 'clim', clim)
    
    hold on;
    plot(cirx, ciry, 'w');

end

figure;
plot(0:M,Err);
xlabel('M');
ylabel('Avg Error');

e = (Err - Err(1))./Err(1)*100;
e2 = (Err(2:end) - Err(1:end-1))./Err(1:end-1)*100;
for m = 2:M+1, 
    e3(m) = min((Err(m:end) - Err(m-1))./Err(m-1)*100); 
end

plot(0:M, e);
hold on;
plot(1:M, e2, 'k');
plot(0:M,e3,'r')

for m = 2:M+1
    if (e2(m) > -0.1)
        Mrec = m - 1;
        break;
    end
end


% for m = 1:M
%     e = (Err(m+1:end) - Err(m))./Err(m)*100;
%     if (mean(e) > -1)
%         Mrec = m - 1;
%         break;
%     end
% end

end

