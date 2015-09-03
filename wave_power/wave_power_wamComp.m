function [aS, aR, comp, waveField, hydBody] = wave_power_wamComp(folder, beam, nTheta, M, varargin)

[opts, args] = checkOptions({{'RndZero'}, {'WaveField', 1}, {'Depth', 1}, ...
    {'HydroBody'}, {'SphereEnd', 1}, {'FlareEnd', 1}, {'CustProf', 1}}, varargin);
rndZero = opts(1);
if (opts(2))
    wfarray = args{2};
else
    wfarray = [];
end

if (opts(3))
    depth = args{3};
else
    depth = 3;
end
compHB = opts(4);

if (opts(5))
    sphLen = args{5};
else
    sphLen = 0;
end

if (opts(6))
    flrRad = args{6};
else
    flrRad = 0;
end

if (opts(7))
    custProf = args{7};
else
    custProf = [];
end

rho = 1000;     
h = 10;          
lam = 10;
T = IWaves.Lam2T(lam, h);

if (compHB)
    N = 40;                         
    beta = 0:2*pi/N:2*pi*(1-1/N);   
else
    beta = 0;
end

% set up flap dimension
rad = 1;
Nr = 6;
Ntheta = 32;
Ny = round(beam/0.17);

% Nr = 2*Nr;
% Ntheta = 2*Ntheta;
% Ny = 2*Ny;

if (sphLen > 0)
    wec = FloatingBristolCyl(rho, rad, beam, depth, Nr, Ntheta, Ny, 'SphereEnd', sphLen);
elseif (flrRad > 0)
    wec = FloatingBristolCyl(rho, rad, beam, depth, Nr, Ntheta, Ny, 'FlareEnd', flrRad);
elseif (~isempty(custProf))
    wec = FloatingBristolCyl(rho, rad, beam, depth, Nr, Ntheta, Ny, 'CustProf', custProf);
else
    wec = FloatingBristolCyl(rho, rad, beam, depth, Nr, Ntheta, Ny);
end
wec.Modes = ModesOfMotion([1 0 1 0 0 0]);   

r = wec.Rcir;                   
dr = 0.1;                       
r = r + dr;                     

nZ = 120;                       

cylArray = BemCylArray(r, nTheta, nZ);                 

wam_run = WamitRunCondition(folder, 'a_run');  

wam_run.Rho = rho;      
wam_run.T = T;                
wam_run.Beta = beta;              
wam_run.H = h;        
wam_run.FloatingBodies = wec;       

wam_run.CylArray = cylArray;   
if (~isempty(wfarray))
    wam_run.FieldArray = wfarray;
    wam_run.ComputeVelocity = true;
end
    
wam_run.WriteRun;          

wam_run.Run;                                    

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;          

waveCir = wam_result.WavePoints;
waveField = wam_result.WaveArray;
hf = wam_result.HydroForces;
comp = HydroBodyComp(hf, wec);

L = 0;
k = (2*pi)./lam;

points = waveCir.FieldPoints;
etaS = waveCir.Elevation('Scattered');
etaR = waveCir.Elevation('Radiated');

[r0, theta, z, eta] = reshapeCirWFPoints(points, etaS{1});
%aS = waveDecomp(k, r0, theta, z, h, eta, L, M, 'SigFigCutoff', 6);
aS = waveDecomp(k, r0, theta, z, h, eta, L, M);

aR = cell(length(etaR), 1);

for n = 1:length(etaR)
    [r0, theta, z, eta] = reshapeCirWFPoints(points, etaR{n});
    aR{n} = waveDecomp(k, r0, theta, z, h, eta, L, M);
end



if (compHB)
    %hydBody = computeHydroBody(waveCir, hf, wec, 'SigFigCutoff', 5,'AccTrim');
    hydBody = computeHydroBody(waveCir, hf, wec, 'SigFigCutoff', 5);
else
    hydBody = [];
end

if (rndZero)
    rndVal = 1e-7;
    
    aMax = max(abs(aS));
    izero = (abs(aS) < rndVal*aMax);
    aS(izero) = 0;
    
    for n = 1:length(aR)
        aMax = max(abs(aR{n}));
        izero = (abs(aR{n}) < rndVal*aMax);
        aR{n}(izero) = 0;
    end
end