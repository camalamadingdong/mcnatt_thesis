function [] = multi_drawBodiesOnPlot(floatBs, varargin)

[opts, args] = checkOptions({'Cir', {'Color', 1}, {'Axis', 1}}, varargin);
drawCir = opts(1);

if (opts(2))
    col = args{2};
else
    col = 'w';
end

if (opts(3))
    ax = args{3};
else
    ax = gca;
end

hold on;

thet = (0:pi/30:2*pi)';
cirth = [cos(thet), sin(thet)];
Nth = length(thet);



for m = 1:length(floatBs)
    cen = floatBs(m).XYpos;
    outline = floatBs(m).WaterPlaneSec;
    Np = size(outline, 1);
        
    ang = floatBs(m).Angle;
    if (ang ~= 0)
        Rot = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];
        for n = 1:Np
            outline(n,:) = Rot*outline(n,:).';
        end
    end
        
    bod = ones(Np, 1)*cen + outline;
    plot(ax, bod(:,1), bod(:,2), col);
    
    if (drawCir)
        R = floatBs(m).Rcir;
        cir = ones(Nth, 1)*cen + R*cirth;
        plot(ax, cir(:,1), cir(:,2), col);
    end
end

end