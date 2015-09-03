function [Ameas, Acyl, Agrid, Awam, am, amW, X, Y, R2c, R2w, Std, incAmp, motAmp] = cylexp_analyzeData(dataDir, body, f, waveType, varargin)

[opts, args] = checkOptions({{'Orient', 1}, {'M', 1}, {'PlotChan', 1}, {'Wam4Ref'}, {'RemoveRef'}, {'IgnoreRuns', 1}}, varargin);
if (opts(1))
    orient = args{1};
else
    orient = 0;
end

if (opts(2))
    M = args{2};
else
    M = 5;
end

if (opts(3))
    plotChans = args{3};
else
    plotChans = {};
end

wam4Ref = opts(4);
removeRef = opts(5);

if (opts(6))
    ignoreRuns = args{6};
else 
    ignoreRuns = [];
end

switch body
    case 'Atten'
        load exp_cyl_hb;
        hydBody = cylHydBody;
    case 'Flap'
        load exp_flap_hb;
        hydBody = flapHydBody;
    otherwise
        error('Body not recognized');
end

switch f
    case 0.8
        fin = 1;
    case 1
        fin = 2;
    case 1.25
        fin = 3;
    otherwise
        error('f not recognized');
end

motAmp = 0;
incAmp = 0;

switch waveType
    case 'Rad'
        isScat = false;
        isTot = false;
                
        [data, wav, amW, motAmp] = computeRadData(dataDir, body, f, plotChans, ignoreRuns, hydBody, M);
        
        xi = data.Motion/180*pi;
        amW = 100*xi*amW;   %scale by body motions and cm..
    case 'Scat'
        isScat = true;
        isTot = false;
        
        [data, datai, wav, amW, ~, incAmp] = computeScatData(dataDir, body, f, orient, plotChans, ignoreRuns, hydBody, M);
    case 'Tot'
        isScat = false;
        isTot = true;
        
        [datar, wavR, amWR] = computeRadData(dataDir, body, f, plotChans, ignoreRuns, hydBody, M);
        [~, datai, wavS, amWS, ami] = computeScatData(dataDir, body, f, orient, plotChans, ignoreRuns, hydBody, M);
        
        data = TestData(dataDir, body, f, 'Motion', 'Free', 'Wave', 'Plane', 'Orient', orient, 'PlotChan', plotChans, 'IgnoreRuns', ignoreRuns);
        
        xi = data.Motion;
        xir = datar.Motion;
        
        amW = amWS + xi/xir*amWR;
        
    otherwise
        error('Wave analysis type not recognized');
end

% Fit coefficients
if (isTot)
    am = wavS.Coefs{1} + xi/xir*wavR.Coefs{1};
    am = {am};
else
    am = wav.Coefs;
end

% Wamit wave
if (isTot)
    wavWR = CirWaves('Out', [0 0], {amWR}, wavR.T, wavR.H);
    wavWS = CirWaves('Out', [0 0], {amWS}, wavS.T, wavS.H);
else
    wavW = CirWaves('Out', [0 0], {amW}, wav.T, wav.H);
end

rho = 1000;

Ameas = cell(6,1);
Acyl = cell(6,1);
Awam = cell(6,1);
Std = cell(6,1);

% Circular points
posc = zeros(24,3);
posc(:,1:2) = data.WGPos('c');

if (isScat)
    Ameasd = data.WaveAmps('c');
    Ameasi = datai.WaveAmps('c');
    Ameas{6} = Ameasd - Ameasi;
    
    Stdd = data.WaveAmps('c', 'Std');
    Stdi = datai.WaveAmps('c', 'Std');
    Std{6} = sqrt((Stdd.^2 + Stdi.^2)./2);
else
    Ameas{6} = data.WaveAmps('c');
    Std{6} = data.WaveAmps('c', 'Std');
end

if (isTot)
    Ameasi = datai.WaveAmps('c');
    wfposWR = CirWaveField(rho, wavWR, false, posc);
    wfposWS = CirWaveField(rho, wavWS, false, posc);
    
    wfposPer = wfposWS + xi/xir*wfposWR;
    Aw = wfposPer.Elevation;
    Aw = Ameasi + Aw{1};
else
    wfposW = CirWaveField(rho, wavW, false, posc);
    Aw = wfposW.Elevation;
    Aw = Aw{1};
end
Awam{6} = Aw;

% Stuff for reflections
cirref = Ameas{6} - Awam{6};
wavref = TestData.CreateCirWaves('In', cirref, f, M);

% m = -M:M;
% 
% f = data.Frequency;
% h = wav.H;
% k = IWaves.SolveForK(2*pi*f, h);
% r = 0.8;
% kr = k*r;
% 
% Hm = besselh(m,2,kr);
% Jm = besselj(m,kr);
% 
% amI2 = 1./Jm.*(amMHm - amW.*Hm);
% % reflected wave
% wavref = CirWaves('In', [0 0], {amI2}, 1./f, h);

if (wam4Ref)
    if (isTot)
        error('Total wave field not set up with this option');
    end
     wfposi = IncCirWaveField(rho, wavref, false, posc);
     wfpos = wfposW + wfposi;
     amI22 = wavref.Coefs;
     Ac = wfpos.Elevation;
     Ac = Ac{1};
else
    if (isTot)
        wfposR = CirWaveField(rho, wavR, false, posc);
        wfposS = CirWaveField(rho, wavS, false, posc);
        
        wfposPer = wfposS + xi/xir*wfposR;
        Ac = wfposPer.Elevation;
        Ac = Ameasi + Ac{1};
    else
        wfpos = CirWaveField(rho, wav, false, posc);
        Ac = wfpos.Elevation;
        Ac = Ac{1};
    end
end

Acyl{6} = Ac;

% Position along spokes
pos = zeros(8, 3);

for m = 1:5;
    wgNames = ['s' num2str(m)];
    if (isScat)
        Ameasd = data.WaveAmps(wgNames);
        Ameasi = datai.WaveAmps(wgNames);
        
        Ameas{m} = Ameasd - Ameasi;
        
        Stdd = data.WaveAmps(wgNames, 'Std');
        Stdi = datai.WaveAmps(wgNames, 'Std');
        Std{m} = sqrt((Stdd.^2 + Stdi.^2)./2);
    else
        Ameas{m} = data.WaveAmps(wgNames);
        Std{m} = data.WaveAmps(wgNames, 'Std');
    end
    
    pos(:,1:2) = data.WGPos(wgNames);
    
    if (isTot)
        Ameasi = datai.WaveAmps(wgNames);
        wfposWR = CirWaveField(rho, wavWR, false, pos);
        wfposWS = CirWaveField(rho, wavWS, false, pos);
        
        wfposPer = wfposWS + xi/xir*wfposWR;
        Aw = wfposPer.Elevation;
        Aw = Ameasi + Aw{1};
    else
        wfposW = CirWaveField(rho, wavW, false, pos);
        Aw = wfposW.Elevation;
        Aw = Aw{1};
    end
    
    Awam{m} = Aw;
    
    if (wam4Ref)
        wfposi = IncCirWaveField(rho, wavref, false, pos);
        wfpos = wfposW + wfposi;
        amI22 = wavref.Coefs;
        Ac = wfpos.Elevation;
        Ac = Ac{1};
    else
        if (isTot)
            wfposR = CirWaveField(rho, wavR, false, pos);
            wfposS = CirWaveField(rho, wavS, false, pos);
            
            wfposPer = wfposS + xi/xir*wfposR;
            Ac = wfposPer.Elevation;
            Ac = Ameasi + Ac{1};
        else
            wfpos = CirWaveField(rho, wav, false, pos);
            Ac = wfpos.Elevation;
            Ac = Ac{1};
        end
    end

    Acyl{m} = Ac;
end

% remove reflections from the data
if (removeRef)
    k = IWaves.SolveForK(2*pi*f, wav.H);
    r = 0.6:0.2:2;
    A = zeros(length(r), 2);
    for n = 1:length(r)
        A(n, 1) = 1./sqrt(k*r(n))*exp(-1i*k*r(n));
        A(n, 2) = 1./sqrt(k*r(n))*exp(1i*k*r(n));
    end

    for n = 1:5
        eta = Ameas{n}.';
        as = A\eta;
        aout = as(1);
        ain = as(2);
        Ameas{n} = aout./sqrt(k*r).*exp(-1i*k*r);
    end
end

% Compute R2 along each spoke
[R2c, R2w] = computeR2(Ameas, Acyl, Awam);

% Grid wave field
x = -3:0.05:2.4;
y = -8:0.05:4.2;
[X, Y] = meshgrid(x, y);

if (wam4Ref)
    wfgridW = CirWaveField(rho, wavW, true, X, Y);
    wfgridi = IncCirWaveField(rho, wavref, true, X, Y);
    wfgrid = wfgridW + wfgridi;
else
    if (isTot)
        wavI = CirWaves('In', [0 0], ami, wavR.T, wavR.H);
        
        wfgridI = IncCirWaveField(rho, wavI, true, X, Y);
        wfgridR = CirWaveField(rho, wavR, true, X, Y);
        wfgridS = CirWaveField(rho, wavS, true, X, Y);
        
        wfgrid = wfgridI + wfgridS + xi/xir*wfgridR;
    else
        wfgrid = CirWaveField(rho, wav, true, X, Y);
    end
end

tank = wfe_tank_points();
wfgrid.RemoveGeometries(tank, 'In');
Ag = wfgrid.Elevation;
Agrid = Ag{1};

end

function [am] = fourierCir(etac, M)
% Creates circular waves (with the associated coefficients) for
% this run from the circular wave gauge array
theta = 0:2*pi/24:2*pi*23/24;
Ntheta = length(theta);

X = 1/Ntheta*fft(etac);

am = zeros(1, 2*M+1);
am(M+1) =  X(1);
am(M+2:2*M+1) = X(2:M+1);
am(M:-1:1) = X(Ntheta:-1:Ntheta-M+1);

end

function [R2c, R2w] = computeR2(Ameas, Acyl, Awam)
% compute R2
Am = zeros(5*8, 1);
Ac = zeros(5*8, 1);
Aw = zeros(5*8, 1);
for m = 1:5
    i1 = (m-1)*8 + 1;
    i2 = i1 + 7;
    Am(i1:i2) = Ameas{m};
    Ac(i1:i2) = Acyl{m};
    Aw(i1:i2) = Awam{m};
end

mAm = mean(Am);
SStot = sum(abs(Am - mAm*ones(5*8,1)).^2);

SSresC = sum(abs(Ac - Am).^2);
SSresW = sum(abs(Aw - Am).^2);

R2c = 1 - SSresC/SStot;
R2w = 1 - SSresW/SStot;

end

function [data, wav, amW, motAmp] = computeRadData(dataDir, body, f, plotChans, ignoreRuns, hydBody, M)

switch f
    case 0.8
        fin = 1;
    case 1
        fin = 2;
    case 1.25
        fin = 3;
    otherwise
        error('f not recognized');
end

data = TestData(dataDir, body, f, 'Motion', 'Rad', 'PlotChan', plotChans, 'IgnoreRuns', ignoreRuns);

motAmp = data.Motion;

etac = data.WaveAmps('c');
amMHm = fourierCir(etac, M);

wav = data.CirWaves(M);
        
amW = hydBody.RadCoefs{fin};
amW = CirWaveComp.Resize2M(amW, M);
end

function [data, datai, wav, amW, ami, incAmp] = computeScatData(dataDir, body, f, orient, plotChans, ignoreRuns, hydBody, M)

switch f
    case 0.8
        fin = 1;
    case 1
        fin = 2;
    case 1.25
        fin = 3;
    otherwise
        error('f not recognized');
end

data = TestData(dataDir, body, f, 'Motion', 'Fix', 'Wave', 'Plane', 'Orient', orient, 'PlotChan', plotChans, 'IgnoreRuns', ignoreRuns);
datai = TestData(dataDir, 'None', f, 'Wave', 'Plane', 'WGSetup', 2, 'PlotChan', plotChans, 'IgnoreRuns', ignoreRuns);

ia = datai.WaveAmps('all');
incAmp = mean(abs(ia));

cd = data.WaveAmps('c');
ci = datai.WaveAmps('c');
etac = cd - ci;

amMHm = fourierCir(etac, M);

wav = TestData.CreateCirWaves('Out', etac, f, M);

wi = datai.CirWaves(M);
ami = wi.Coefs;
if (orient ~= 0)
    hydBody.Angle = orient;
end
DTM = hydBody.DiffTM{fin};
DTM = CirWaveComp.Resize2M(DTM, M);

amW = DTM*ami{1}.';
amW = amW.';

end
