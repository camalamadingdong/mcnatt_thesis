%% Folder

folder = [myMatPath '\thesis\cylexp'];  
rho = 1000;    
h = 1.16;      
f = [0.8, 1, 1.25];
T = 1./f;

%% Make Cylinder

run_name = 'exp_cyl'; 

% set up cylinder dimension
len = 0.8;
rad = 0.08;
sphereRad = 0.08;
hinPos = [];
Nx = 100;
Ntheta = 32;

wec = FloatingSphereEndCyl(rho, len, rad, sphereRad, hinPos, Nx, Ntheta);
wec.Modes = ModesOfMotion([0 0 0 0 1 0]);   

%% Cylinder Hydrobody        

N = 32;  
Beta = 0:2*pi/N:2*pi*(1-1/N);   

r = wec.Rcir;                   
dr = 0.1;                       
r = r + dr;                     

nZ = 101;                       
nTheta = 2^7; 

cylArray = BemCylArray(r, nTheta, nZ);

wam_run = WamitRunCondition([folder '\wam_runs'], run_name);  

wam_run.Rho = rho;      
wam_run.T = T;                
wam_run.Beta = Beta;              
wam_run.H = h;        
wam_run.FloatingBodies = wec;       

wam_run.CylArray = cylArray;
wam_run.WriteRun;               
                       
wam_run.Run;           

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;          

waveCir = wam_result.WavePoints;
hydroForces = wam_result.HydroForces;

cylHydBody = computeHydroBody(waveCir, hydroForces, wec, 'SigFigCutoff', 5,...
    'AccTrim');

save([folder '\cylexp_cyl_hb'], 'cylHydBody');

%% Cylinder Grid
                     
Beta = [0 7/4*pi];    

gridPnts = BemFieldArray([-2 -2 0], [0.04 0.04 1], [101 101 1]);

wam_run = WamitRunCondition([folder '\wam_runs'], run_name);  

wam_run.Rho = rho;      
wam_run.T = T;                
wam_run.Beta = Beta;              
wam_run.H = h;        
wam_run.FloatingBodies = wec;       

wam_run.FieldArray = gridPnts;
wam_run.WriteRun;               
                       
wam_run.Run;           

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;          

cylWF = wam_result.WaveArray;
hydroForces = wam_result.HydroForces;

cylComp = HydroBodyComp(hydroForces, wec);

save([folder '\cylexp_cyl_grid'], 'cylComp', 'cylWF');

%% Make Flap

run_name = 'exp_flap';         

len = 0.08;
beam = 0.6;
draft = 0.4;
Nx = 5;
Ny = 38;
Nz = 25;

wec = FloatingFlap(rho, len, beam, draft, Nx, Ny, Nz);
wec.Modes = ModesOfMotion([0 0 0 0 1 0]);   

%% Flap Hydrobody

N = 32;                         
Beta = 0:2*pi/N:2*pi*(1-1/N);   

r = wec.Rcir;                   
dr = 0.1;                       
r = r + dr;                     

nZ = 101;                       
nTheta = 2^7; 

cylArray = BemCylArray(r, nTheta, nZ);

wam_run = WamitRunCondition([folder '\wam_runs'], run_name);  

wam_run.Rho = rho;      
wam_run.T = T;                
wam_run.Beta = Beta;              
wam_run.H = h;        
wam_run.FloatingBodies = wec;       

wam_run.CylArray = cylArray;
wam_run.WriteRun;              

wam_run.Run; 

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;          

waveCir = wam_result.WavePoints;
hydroForces = wam_result.HydroForces;

flapHydBody = computeHydroBody(waveCir, hydroForces, wec, 'SigFigCutoff', 5,...
    'AccTrim');

save([folder '\cylexp_flap_hb'], 'flapHydBody');

%% Flap Grid
                     
Beta = [0 7/4*pi];    

gridPnts = BemFieldArray([-2 -2 0], [0.04 0.04 1], [101 101 1]);

wam_run = WamitRunCondition([folder '\wam_runs'], run_name);  

wam_run.Rho = rho;      
wam_run.T = T;                
wam_run.Beta = Beta;              
wam_run.H = h;        
wam_run.FloatingBodies = wec;       

wam_run.FieldArray = gridPnts;
wam_run.WriteRun;               
                       
wam_run.Run;           

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;          

flapWF = wam_result.WaveArray;
hydroForces = wam_result.HydroForces;

flapComp = HydroBodyComp(hydroForces, wec);

save([folder '\cylexp_flap_grid'], 'flapComp', 'flapWF');

%% Flap Grid