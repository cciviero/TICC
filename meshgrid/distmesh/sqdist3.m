function [d] = sqdist3(xx,yy,zz,xs,ys,zs)
% function [d] = sqdist3(xx,yy,zz,xs,ys,zs)
  
d = (xx-xs).^2 + (yy-ys).^2 + (zz-zs).^2; 
