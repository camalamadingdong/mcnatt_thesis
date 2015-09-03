function [Aa, Ba, Fexa, Xia] = multi_computeHBdist(hydBod, Dist, beta, side2side, skipDist)

hb1 = HydroBody(hydBod);
hb1.XYpos = [0 -100];
hb2 = HydroBody(hydBod);
hb2.XYpos = [0 100];

for n = 1:length(beta)
    iwaves(n) = PlaneWaves(ones(size(hydBod.T)), hydBod.T, beta(n)*ones(size(hydBod.T)), hydBod.H);
end

hcomp = HydroArrayComp([hb1 hb2], iwaves);
%hcomp = HydroArrayComp_jfmScale([hb1 hb2]);
dof = hcomp.DoF;
nT = length(hcomp.T);

nbeta = length(beta);

ndist = length(Dist);

Aa = zeros(ndist, nT, dof, dof);
Ba = zeros(ndist, nT, dof, dof);
Fexa = zeros(ndist, nT, nbeta, dof);
Xia = zeros(ndist, nT, nbeta, dof);

for n = 1:ndist
    dist = Dist(n);
    if (dist <= skipDist)
        Aa(n,:,:,:) = NaN*ones(nT, dof, dof);
        Ba(n,:,:,:) = NaN*ones(nT, dof, dof);
        Fexa(n,:,:,:) = NaN*ones(nT, nbeta, dof);
        Xia(n,:,:,:) =  NaN*ones(nT, nbeta, dof);
    else
        if (side2side)
            hcomp.BodXY = [0 -dist/2; 0 dist/2];
        else
            hcomp.BodXY = [-dist/2 0; dist/2 0];
        end
        Aa(n,:,:,:) = hcomp.A;
        Ba(n,:,:,:) = hcomp.B;
        Fexa(n,:,:,:) = hcomp.Fex;
        Xia(n,:,:,:) =  hcomp.Motions;
    end
end