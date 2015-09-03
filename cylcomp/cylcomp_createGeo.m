function [wec] = cylcomp_createGeo(d, type, varargin)
    

rho = 1000;
g = IWaves.G;

switch type
    case 'Cyl'
        draft = 3/2*d;  

        Ntheta = 48;
        Nr = 8;
        Nz = 24;

        wec = FloatingCylinder(rho, d/2, draft, draft, Ntheta, Nr, Nz);

        wec.Handle = 'cyl';                                            

        % Dh = 180
        % Kh = 1100
        D = zeros(6, 6);
        K = zeros(6, 6);
        dimD = rho*sqrt(g)*d^5/2;
        dimK = rho*g*d^2;

        if (~isempty(varargin))
            motion = varargin{1};
            switch motion
                case 'Surge'
                    wec.Modes = ModesOfMotion([1 0 0 0 0 0]);
                    wec.Handle = 'surge';
                    Dnd = 1.852;
                    Knd = 1.071;
                    D(1,1) = dimD*Dnd;
                    K(1,1) = dimK*Knd;
                case 'Heave'
                    wec.Modes = ModesOfMotion([0 0 1 0 0 0]); 
                    wec.Handle = 'heave';
                    Dnd = 0.115;
                    Knd = 0.112;
                    D(3,3) = dimD*Dnd;
                    K(3,3) = dimK*Knd;
                otherwise
                    error('Cylinder type not recognized');
            end
        end
    case 'Flap'
        len = d/4;
        beam = 4*d;
        draft = (3*pi/8)*d;

        Nx = 4;
        Ny = 60;
        Nz = 18;

        wec = FloatingFlap(rho, len, beam, draft, Nx, Ny, Nz, 'NoInt');
        wec.Modes = ModesOfMotion([0 0 0 0 1 0]);
        wec.Handle = 'flap';                                   

        D = zeros(6, 6);
        K = zeros(6, 6);
        dimD = rho*sqrt(g)*d^9/2;
        dimK = rho*g*d^4;
        Dnd = 0.607;
        Knd = 5.303;
        D(5,5) = dimD*Dnd;
        K(5,5) = dimK*Knd;

        wec.Dpto = D;
        wec.K = K;
    case 'Atten'
        conePct = 0.05;

        len = 10*d;
        %len = 8*d;
        %beam = sqrt(45/104)*d;
        beam = sqrt(45/13/len);
        radius = beam/2;

        Nx = 120;
        Ntheta = 16;

        wec = FloatingAttenuator(rho, len, radius, conePct, Nx, Ntheta, 'NoInt');   
        wec.Handle = 'atten';                                            

        wec.Modes.Vector = [0 0 1 0 0 0 1];

        D = zeros(7, 7);
        K = zeros(7, 7);
        % d = 19100
        dimD = rho*sqrt(g)*d^9/2;
        dimK = rho*g*d^4;
        Dnd = 12.2;
        Knd = 0;
        D(7,7) = dimD*Dnd;
        K(7,7) = dimK*Knd;
    otherwise
        error('Geometry type not recognized');
end

wec.Dpto = D;
wec.K = K;

end