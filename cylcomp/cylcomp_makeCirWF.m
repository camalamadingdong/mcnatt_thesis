function [cWave, fWave, r0, Bs, fFuncs] = cylcomp_makeCirWF(waveC, N, M, X, Y)

points = waveC.FieldPoints;

h = waveC.H;
rho = waveC.Rho;
T = waveC.T;
beta = waveC.IncWaveVals;

% Scattered
eta = waveC.Elevation('Scattered');
[r0, theta, z, eta] = reshapeCirWFPoints(points, eta{1});    
ks = IWaves.SolveForK(2*pi/T, h, 'Evanescent', N);
B = waveDecomp(ks, r0, theta, z, h, eta, N, M, 'SigFigCutoff', 5);

cwave = CirWaves('Out', [0 0], {B}, T, h);
cWfS = CirWaveField(rho, cwave, true, X, Y, 'NoVel');
cWfS = WaveFieldCollection(cWfS, 'Direction', beta);

fWfS = KochinWaveField(rho, cwave, true, X, Y, 'NoVel');
fWfS = WaveFieldCollection(fWfS, 'Direction', beta);

pwave = PlaneWaves(1, T, beta, h);
wfI = PlaneWaveField(rho, pwave, true, X, Y);
wfI = WaveFieldCollection(wfI, 'Direction', beta);

% Radiated
etaR = waveC.Elevation('Radiated');
nR = length(etaR);

Bs = cell(1, 1 + nR);
Bs{1} = B;

fFuncs = cell(1, 1 + nR);
fFuncs{1} = cwave.KochinFunc;

for n = 1:nR
    [r0, theta, z, eta] = reshapeCirWFPoints(points, etaR{n});
    ks = IWaves.SolveForK(2*pi/T, h, 'Evanescent', N);
    B = waveDecomp(ks, r0, theta, z, h, eta, N, M, 'SigFigCutoff', 5);
    
    cwave = CirWaves('Out', [0 0], {B}, T, h);
    
    cWfR = CirWaveField(rho, cwave, true, X, Y, 'NoVel');
    cWfRs(n) = cWfR; 
    
    fWfR = KochinWaveField(rho, cwave, true, X, Y, 'NoVel');
    fWfRs(n) = fWfR;
    
    Bs{1+n} = B;
    fFuncs{1+n} = cwave.KochinFunc;
end

cWfR = WaveFieldCollection(cWfRs);
fWfR = WaveFieldCollection(fWfRs);

cWave = FBWaveField(wfI, cWfS, cWfR);
fWave = FBWaveField(wfI, fWfS, fWfR);
