function [comp, wec, hydBod, waveField] = multi_makeRun(wecType, runType, runFold, lam, beta, varargin)

[opts, args] = checkOptions({{'Scale',1}, {'Front2Back'}, {'Side2Side'}, {'Symmetry'}, {'Rotate', 1}}, varargin);
if(opts(1))
    scale = args{1};
else
    scale = 1;
end

if (opts(2))
    frontBack = true;
elseif (opts(3))
    frontBack = false;
else
    frontBack = false;
end

if (opts(4))
    useSym = true;
else
    useSym = false;
end

if (opts(5))
    rotate = true;
    angs = args{5};
    izero = [];
else
    rotate = false;
end

if (strcmp(runType, 'Dist'))
    dist = varargin{1};
end

if (strcmp(runType, 'WaveField'))
    fieldArray = varargin{1};
end

if (strcmp(runType, 'MedArray'))
    pos = varargin{1};
    fieldArray = varargin{2};
end

h = 10;    
rho = 1000;


switch wecType
    case 'Cyl'
        hei = scale*1.5;
        draft = scale*1;

        Ntheta = 64;
        Nr = 8;
        Nz = 12;

        wec = FloatingCylinder(rho, 1, hei, draft, Ntheta, Nr, Nz, 'NoInt');

        wec.Handle = 'cyl';     
        wec.Modes = ModesOfMotion([1 1 1 1 1 1]);
    case 'Atten'
        len = scale*10;
        rad = scale*0.5;

        Nx = 120;
        Ntheta = 16;

        sphereRad =0.1*len;
        hingePos = [-len/6; len/6];

        wec = FloatingSphereEndCyl(rho, len, rad, sphereRad, hingePos, Nx, Ntheta, 'Notch');  
        wec.Handle = 'atten';  
    otherwise
        wec = [];
end

comp = [];
hydBod = [];

if (~strcmp(runType, 'Geo'))
    name = 'a_run';

    wam_run = WamitRunCondition(runFold, name); 
    wam_run.MaxItt = 100;
    wam_run.Rho = rho;   
    
    wam_run.FloatingBodies = wec;       

    T = IWaves.Lam2T(scale*lam, h);
    wam_run.T = T;
    
    switch runType
        case 'HydroBody'
            N = 80;

            if (useSym)
                switch wecType
                    case 'Cyl'
                        beta = 0;
                    case 'Atten'
                        beta = 0:2*pi/N:pi;
                end
            else
                beta = 0:2*pi/N:2*pi*(1-1/N);
            end

            dr = 0.1;
            nZ = 201;
            nTheta = 2^8;

            r = wec.Rcir + dr;

            wam_run.CylArray = BemCylArray(r, nTheta, nZ);
        case 'Dist'
            wec2 = FloatingBody(wec);
        
            if (frontBack)
                wec.XYpos = [-dist/2 0];
                wec2.XYpos = [dist/2 0];
            else
                wec.XYpos = [0 -dist/2];
                wec2.XYpos = [0 dist/2];
            end

            wam_run.FloatingBodies = [wec wec2];
        case 'WaveField'
            wam_run.FieldArray = fieldArray;
        case 'MedArray'
            Nwec = size(pos,1);
            for m = 1:Nwec
                wecs(m) = FloatingBody(wec);
                wecs(m).Handle = ['Wec' num2str(m)];
                wecs(m).XYpos = pos(m,:);
                if (rotate)
                    wecs(m).Angle = angs(m);
                    if (isempty(izero))
                        if (angs(m) == 0)
                            izero = m;
                        end
                    end
                        
                end
            end
            wam_run.FloatingBodies = wecs;       
            
            wam_run.FieldArray = fieldArray; 
            
            wam_run.NCPU = 1;
            wam_run.RAMGBmax = 5;
    end
    wam_run.Beta = beta;                                   
    wam_run.H = h;

    % wam_run.NCPU = 3;
    % wam_run.RAMGBmax = 5;

    wam_run.WriteRun;   

    wam_run.Run; 

    wam_result = WamitResult(wam_run);                   
    wam_result.ReadResult;                              

    hydroF = wam_result.HydroForces;
    floatB = wam_result.FloatingBodies;
    waveField = wam_result.WaveArray;
    
    if (rotate)
        % due to a bug the hydrostatic force out of WAMIT is not correct
        % for rotated bodies in an array - need to do another computation
        % for an unrotated body
        
         name = 'a_run';

        wam_run = WamitRunCondition(runFold, name); 
        wam_run.MaxItt = 100;
        wam_run.Rho = rho;   

        wam_run.FloatingBodies = wec;       
        wam_run.T = 1;
            
        wam_run.NCPU = 1;
        wam_run.RAMGBmax = 5;
        wam_run.Beta = 0;                                   
        wam_run.H = h;
        
        wam_run.WriteRun;   
        wam_run.Run; 
        
        wam_result = WamitResult(wam_run);                   
        wam_result.ReadResult;  
        
        C0 = wam_result.HydroForces.C;
        
        dof = wec.Modes.DoF;
        C = zeros(Nwec*dof, Nwec*dof);

        for n = 1:Nwec
            ind1 = (n-1)*dof+1;
            ind2 = ind1+dof-1;
            C(ind1:ind2, ind1:ind2) = C0;
        end

        hydroF = HydroForces(hydroF.T, hydroF.Beta, hydroF.A, hydroF.B, C, hydroF.Fex, hydroF.H, hydroF.Rho);
    end
    
    switch wecType
        case 'Cyl'
            comp = HydroBodyComp(hydroF, floatB);
        case 'Atten'
            %comp = HydroBodyComp(hydroF, floatB);
            comp = HydroBodyComp(hydroF, floatB, 'UseBodC');
    end
    
    if (strcmp(runType, 'HydroBody'))
        waveC = wam_result.WavePoints;
        
        if (useSym)
            switch wecType
                case 'Cyl'
                    hydBod = computeHydroBody(waveC, hydroF, floatB, 'SigFigCutoff', 5, 'AxisSym', wec.Modes);
                case 'Atten'
                    hydBod = computeHydroBody(waveC, hydroF, floatB, 'SigFigCutoff', 5, 'XSym');
            end
        else
            hydBod = computeHydroBody(waveC, hydroF, floatB, 'SigFigCutoff', 5, 'AccTrim');
        end
    end
end

