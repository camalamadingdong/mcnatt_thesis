function [] = cylexp_plotWg(wgSetup, varargin)

ctheta = linspace(0,2*pi,101);
r = 0.8;

hold on;
plot(r*cos(ctheta), r*sin(ctheta), varargin{:});

if (wgSetup == 1)
    plot([-2*r 2*r], [0 0], varargin{:});
    for n = 5:7
        plot(2*r*[cos(n*pi/4), 0], 2*r*[sin(n*pi/4), 0], varargin{:});
    end
elseif (wgSetup == 2)
    for n = 4:8
        plot(cos(n*pi/4)*[0.6, 2], sin(n*pi/4)*[0.6, 2], varargin{:});
    end
else
end

end