% =========================================================================
%  FILENAME
%          marias_postproc.m
%  AUTHOR
%          just adjust Karins code by Maria for RHUM-RUM
%  PURPOSE
%          distmesh: Global tomography mesh
%  INFORMATION
%          
%  INPUT
%        
%  OUTPUT 
%         
%  NOTE   
%        
% =========================================================================

%% ------------------ PART 1: SET PARAMETERS

meshname = 'RR1115_v2'; 
hmatrix_mode = 'stationdensified';
fstatlist = 'station_file.txt';
sig = 60;

dim=3;  % dimensions  

rplus = 50; 
R0 = 6371 + rplus;

%% ------------------ PART 2: PATHS AND NAMES
basedir = ['/Users/maria/PhD/Codes/TICC/meshgrid'];

fout_facets   = sprintf('facets.%s',meshname);
fout_vertices = sprintf('vertices.%s',meshname);
fout_resl     = sprintf('mesh.resolvinglength.%s',meshname);

% CD to workdir
workdir = fullfile(basedir,meshname);
if(not(exist(workdir,'dir')))
  [success,msg,msgid] = mkdir(workdir);
  if(not(success)) error(sprintf('%s %s',msg,msgid)); end
end  
cd(workdir);

% Path to distmesh software
distmeshpath = ['/Users/maria/PhD/Codes/TICC/meshgrid/distmesh'];

addpath(distmeshpath);

%% ------------------ PART 3: FUNCTIONS

%--------------------------------------------
% Distance function: Earth's shape
%--------------------------------------------
expr = sprintf('sqrt(sum(p.^2,2))-%f',R0);
d = inline(expr,'p');
pfix = [];  % fixed points

%--------------------------------------------
% Resolution width function hh 
%--------------------------------------------
% use matrix h to interpolate
delta = 100;  % 
RX = 1.1*R0;
[xx,yy,zz] = meshgrid(-RX:delta:RX,-RX:delta:RX,-RX:delta:RX);
depth = 6371 - sqrt(sqdist3(xx,yy,zz,0,0,0));

hmax = 600;  % largest edge length ever

% define depth dependence of h as piecewise
% linear splines between d1, d2, and d3.
% this part is important because it defines how the mesh is calculated
d1 = 0;     h1 = 400;
d2 = 1000;  h2 = 400; % regular grid in the upper 1000km 
d3 = 2889;  h3 = hmax; % gird is getting linear bigger towards CMB

a12 = (h1-h2)/(d1-d2); b12 = h1-a12*d1;
a23 = (h2-h3)/(d2-d3); b23 = h2-a23*d2;

hh = hmax*ones(size(zz));

ii = find(depth>-100 & depth<=d2);
hh(ii)  = a12*depth(ii)+b12;

ii = find(depth>d2 & depth<=d3);
hh(ii)  = a23*depth(ii)+b23;


%% ------------------ PART 4: POST-PROCESSING


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
ucall = sprintf('mv %s %s.bkp; qhull QJ d i < %s > %s',fnam2,fnam2,fnam1,fnam2);
unix(ucall);

if(0)
% OUTFILE: facets.makemymesh
% Write facets file consistent with the way qhull.f
% expects it, e.g. edge numbering starts at 0 
% not 1 (like Matlab's delaunayn). Just subtract 1
% from each t
fnam = 'facets.makemymesh';
disp(sprintf('Writing facets file: %s\n',fnam));
fid = fopen(fnam,'w');
fprintf(fid,'%12d\n',length(t));
for kk=1:length(t)
  fprintf(fid,'%12d %12d %12d %12d\n',...
	    t(kk,1)-1,t(kk,2)-1,t(kk,3)-1,t(kk,4)-1);
end  
fclose(fid);
end % if(0) 

%-----------------------

% Resolution length at every vertex p
% fh = 'hmatrix3d';  
% HH0=feval(fh,p(:,:),xx_p,yy_p,zz_p,hh,hh);

% OUTFILE: mesh.resolvinglength (format like solx.* files)
% First line is a dummy line
%  3 121.59749 607.98743   0.20000   0.30000 dlnVp 39932
%  write(9,fmt='(i3,4f10.5,1x,a5)') 3,epsnorm,epssmooth,epsratio,
%     &  ytry(itry),chcor5

% fnam = fout_resl;
% % fnam = 'mesh.resolvinglength';
% fid = fopen(fnam,'w');
% fprintf(fid,'%3d%10.5f%10.5f%10.5f%10.5f %5s\n',...
% 	3,100,500,0.2,0.3,'dlnVp');  % bogus line, simulates solx format
% fprintf(fid,'%12d\n',length(p));
%   for kk=1:length(p)
%     fprintf(fid,'%12.4f %12.4f %12.4f %12.4f\n',...
% 	    p(kk,1),p(kk,2),p(kk,3),HH0(kk));
%   end  
% fclose(fid);


%-------   PLOTS  -----------------

% Viz results: half-sphere at different angles
% alpha = tan(phi) = y/x = p(:,2)/p(:,1)

phicut = 90;  % cut/view longitude in degrees
expr = sprintf('p(:,2)- %f*p(:,1)<0',tan(phicut)*(pi/180));
expr = 'p(:,1)>0';
expr = 'p(:,1)>0';

fignam = sprintf('fig_meshcut.%s',meshname)
hfig = 555; figure(hfig);
simpplot(p,t,expr);
% simpplot(p,t,'p(:,2)>0');
title(sprintf('Final mesh, cut at longitude %d degrees',phicut));
drawnow
saveas(hfig,fignam,'m');


% Histogram of tetrahedra quality
fignam = sprintf('fig_histo.%s',meshname);
hfig = 9999; figure(hfig);
Q=simpqual(p,t,1);
hist(Q,40);
title('Histogram: simplex quality');
saveas(hfig,fignam,'m');



%---  tmp tmp tmp!!  -----------------
% was hacked for USA2
%% basedir = '/raid/AmplitudeProjects/testdistmesh/meshresolution/';
%% unix(sprintf('cp mesh.resolvinglength %s',basedir));
%% unix(sprintf('cp facets.makemymesh %s',basedir));
%% unix(sprintf('cp vertices.makemymesh %s',basedir));

%---  tmp tmp tmp!!  -----------------




%===============   2-D   ================

if(0)

dim=2;   

R0 = 6421;  % 6371

% distance function: Earth's shape
% d=inline('sqrt(sum(p.^2,2))-6421','p');
d=inline('sqrt(sum(p.^2,2))-6421','p');
pfix = [0 0];

% element width function
% interpolate from matrix
delta = 50;
[xx,yy] = meshgrid(-1.1*R0:delta:1.1*R0,-1.1*R0:delta:1.1*R0);
depth = 6371 - sqrt(xx.^2+yy.^2);
[i1, i2] = size(depth);
depth = reshape(depth,i1*i2,1);

% define depth dependence of h as piecewise
% linear splines.
hmax = 600;

d1 = 0;    h1 = 100;
d2 = 670;  h2 = 200;
d3 = 2889; h3 = 600;

a12 = (h1-h2)/(d1-d2); b12 = h1-a12*d1;
a23 = (h2-h3)/(d2-d3); b23 = h2-a23*d2;

hh = hmax*ones(size(yy));

ii = find(depth>-100 & depth<=d2);
hh(ii)  = a12*depth(ii)+b12;

ii = find(depth>d2 & depth<=d3);
hh(ii)  = a23*depth(ii)+b23;

hh = reshape(hh,[i1 i2]);


%--------------------------------------

h0 = 100;  % h0, sets repulsive forces, set to average of h(x,y,z) approximately?


[p,t]=distmeshnd_hvar(d,@hmatrix,h0,[-R0*ones(1,dim);R0*ones(1,dim)],pfix,xx,yy,hh,hh);
 
 
end % if(0)
