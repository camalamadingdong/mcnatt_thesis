cmap = colormap(gca);
ncol = size(cmap,1);
cstep = floor(ncol/10);
imap = 1:cstep:ncol;
cmap2 = cmap(imap,:);
colormap(cmap2)