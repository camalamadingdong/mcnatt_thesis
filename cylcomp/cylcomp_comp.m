function [comp, wf1, wf2] = cylcomp_comp(folder, d, lam, beta, type, wftype, varargin)
        
h = 10*d;     
rho = 1000;
name = 'thecomp';

switch type
    case 'Heave'
        wec = cylcomp_createGeo(d, 'Cyl', 'Heave');
    case 'Surge'
        wec = cylcomp_createGeo(d, 'Cyl', 'Surge');
    case 'Flap'
        wec = cylcomp_createGeo(d, 'Flap');
    case 'Atten'
        wec = cylcomp_createGeo(d, 'Atten');
    case 'Atten2'
        wec = cylcomp_createGeo(d, 'Atten2');
    otherwise
        error('Comp type not recognized');
end

wam_run = WamitRunCondition(folder, name);           
wam_run.Rho = rho;                   

wam_run.FloatingBodies = wec;               

lam = d*lam;
wam_run.T = IWaves.Lam2T(lam, h);
wam_run.Beta = beta;                                   
wam_run.H = h;


switch wftype
    case 'Grid'
        gwf = true;
        cwf = false;
    case 'Circ'
        gwf = false;
        cwf = true;
    case 'Both'
        gwf = true;
        cwf = true;
    otherwise
        gwf = false;
        cwf = false;
end

if (gwf)
    wfdist = varargin{1};
    wfdist = d*wfdist;
    del = wfdist/100; 
    
    fieldArray = BemFieldArray([-wfdist -wfdist 0], [del del 1], [201 201 1]);
    wam_run.FieldArray = fieldArray;
end
 
if (cwf)
    r = wec.Rcir;         
    dr = 0.1*d;
    r = r + dr;                     

    nZ = 200;                       
    nTheta = 2^8;  
    cylArray = BemCylArray(r, nTheta, nZ);         
    wam_run.CylArray = cylArray; 
end

wam_run.WriteRun;    

wam_run.Run;

wam_result = WamitResult(wam_run);                 
wam_result.ReadResult;                              

hydroForces = wam_result.HydroForces;
wec = wam_result.FloatingBodies;                  
comp = HydroBodyComp(hydroForces, wec);  
                 
switch wftype
    case 'Grid'
        wf1 = wam_result.WaveArray;  
        wf2 = [];
    case 'Circ'
        wf1 = wam_result.WavePoints;
        wf2 = [];
    case 'Both'
        wf1 = wam_result.WaveArray;  
        wf2 = wam_result.WavePoints;
    otherwise
        wf1 = [];
        wf2 = [];
end


end

