% Reload mesh (p,t) if run was interrupted
load dmeshtmp.mat  % recover p t (if run wasn;t fully converged)
% OUTFILE: vertices.makemymesh
fnam1 = fout_vertices;
% fnam1 = 'vertices.makemymesh';
disp(sprintf('Writing vertices file: %s\n',fnam1));
fid = fopen(fnam1,'w');
fprintf(fid,' 3\n');
fprintf(fid,'%12d\n',length(p));
for kk=1:length(p)
fprintf(fid,'%12.4f %12.4f %12.4f\n',...
p(kk,1),p(kk,2),p(kk,3));
end
fclose(fid);


% OUTFILE: facets.makemymesh
fnam2 = fout_facets;
% fnam2 = 'facetscd me.makemymesh';
% ucall = sprintf('mv %s %s.bkp; qhull QJ d i < %s > %s',fnam2,fnam2,fnam1,fnam2);
ucall = sprintf('qhull QJ d i < %s > %s', fnam1, fnam2);
unix(ucall);

fnam = 'facets.makemymesh';
disp(sprintf('Writing facets file: %s\n',fnam));
fid = fopen(fnam,'w');
fprintf(fid,'%12d\n',length(t));
for kk=1:length(t)
fprintf(fid,'%12d %12d %12d %12d\n',...
t(kk,1)-1,t(kk,2)-1,t(kk,3)-1,t(kk,4)-1);
end
fclose(fid);

fh = 'hmatrix3d';
HH0=feval(fh,p(:,:),xx,yy,zz,hh,hh);
% OUTFILE: mesh.resolvinglength (format like solx.* files)
% First line is a dummy line
%  3 121.59749 607.98743   0.20000   0.30000 dlnVp 39932
%  write(9,fmt='(i3,4f10.5,1x,a5)') 3,epsnorm,epssmooth,epsratio,
%     &  ytry(itry),chcor5
fnam = fout_resl;
% fnam = 'mesh.resolvinglength';
fid = fopen(fnam,'w');
fprintf(fid,'%3d%10.5f%10.5f%10.5f%10.5f %5s\n',...
3,100,500,0.2,0.3,'dlnVp');  % bogus line, simulates solx format
fprintf(fid,'%12d\n',length(p));
for kk=1:length(p)
fprintf(fid,'%12.4f %12.4f %12.4f %12.4f\n',...
p(kk,1),p(kk,2),p(kk,3),HH0(kk));
end
fclose(fid);


phicut = 0;  % cut/view longitude in degrees
expr = sprintf('p(:,2)- %f*p(:,1)<0',tan(phicut)*(pi/180));
expr = 'p(:,1)>0';
fignam = sprintf('fig_meshcut.%s',meshname)
hfig = 555; figure(hfig);
simpplot(p,t,expr);
% simpplot(p,t,'p(:,2)>0');
title(sprintf('Final mesh, cut at longitude %d degrees',phicut));
drawnow

