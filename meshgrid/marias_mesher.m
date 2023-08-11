% =========================================================================
%  FILENAME
%          marias_mesher.m
%  AUTHOR
%          just adjust Karins code by Maria for RHUM-RUM
%  PURPOSE
%          distmesh: Global tomography mesh
%  INFORMATION
%          It CALLS the following:
%          * functions from package distmesh (by Persson & Strang, SIAM Review,
%          2004), incl. simpplot.m, simpqual.m
%          distmeshnd_hvar.m, written by Karin as generalized version of
%          distmeshnd.m to read to accept 3-D edge-contraint matrix H as 
%          as argument.
%          Aux function sqdist3.m (Karin)
%          * qhull, which must be running under linux
%          * Guust's software 
%          /data/jingjie/sigloch/INVERSION/Programs/buildgrid/sol2xyz
%          for the case hmatrix_mode = 'kerneldensified', which is disabled
%          for the time being. So that software is not needed.
%  INPUT
%          * parameters set in the first section below
%          * the station list text file specified by "fstatlist"
%          (provided that hmatrix_mode = 'stationdensified')
%          * the kernel column density file if 
%
%  OUTPUT 
%         (locally to "workdir"):
%          vertex file
%          facet file
%          work directories in the past included:
%          /raid/AmplitudeProjects/testdistmesh
%          /raid/Private/sigloch/testdistmesh
%  NOTE   
%         it may take extremely long or forever for distmeshnd_hvar to
%         reach its break condition. Instead, it can be interrupted whenever the
%         mesh looks good enough (see plot that gives output every 5 iters). Every
%         5th iter also gets saved: save dmeshtmp p t

% =========================================================================

% maybe some dummy stations are necessary 
%             -90.0000              0.00         MT.DUMY.OO.BHZ             0.0        0.0
%              30.0000              -120.00       MT.DUMY.O1.BHZ             0.0        0.0
%               0.0000              0.00         MT.DUMY.O2.BHZ             0.0        0.0
%             -20.0000              -160.00       MT.DUMY.O3.BHZ             0.0        0.0



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

%--------

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


%% ------------------ PART 4: MAIN PART

%-------------------------------------
% Influence of individual stations on
% local resolving length
%-------------------------------------

switch hmatrix_mode 
  %-----------------------
  case 'stationdensified'
  %-----------------------

    % decrease targeted edge length h in the vicinity of 
    % every station read from list. vicinity = a Gaussian
    % of radius sig (in km).
    % Computes matrix h (values [0,1]) and divides
    % hh, the density matrix obtained above:
    % hh = hh./h;
    
    fnam = fstatlist;
    [stla,stlo]=textread(fnam,'%f%f%*[^\n]');
%     [stla,stlo]=unique([stla stlo],'rows');
    nstat = length(stla);

    theta = (90-stla(:))*(pi/180);
    phi   = stlo(:)*(pi/180); 
    rho   = 6371;

    zs    = rho.*cos(theta);
    xs    = rho*sin(theta).*sin(phi);
    ys    = rho*sin(theta).*cos(phi);
    % (Note: the expected expression below messes up the phase)
    %  xs = rho*sin(theta).*cos(phi);  % orig, but problem
    %  ys = rho*sin(theta).*sin(phi);  % orig, but problem
    %    MESHGRID is like NDGRID except that the order of the first two input
    %    and output arguments are switched (i.e., [X,Y,Z] = MESHGRID(x,y,z)
    %    produces the same result as [Y,X,Z] = NDGRID(y,x,z)).  Because of
    %    this, MESHGRID is better suited to problems in cartesian space,
    %    while NDGRID is better suited to N-D problems that aren't spatially
    %    based.  MESHGRID is also limited to 2-D or 3-D.
    %% but interpn (as used in hmatrix3d) expects that the grid was
    %% generated with NDGRID, not meshgrid.
    %% 

    %----------------------------------------------
    % CHOOSE A WEIGHTING FUNCTION around stations
    % and its parameters
    %----------------------------------------------

    % Linear weighting produces more even, smoothly varying 
    % edge lengths. Was used for grids USA5 and later.
    
    
    linear_weighting = 1;
    if(linear_weighting)
      % 66 km min facet length: used in USA10, G10
      hval0 =  66;  % in km (target edge length at 0 km radius, i.e. the station)
      % 100 km min facet length: used in G11
      % hval0 =  100;  % in km (target edge length at 0 km radius, i.e. the station)

      R = 1200;  % km radius
      hvalR = 200;  % edge length in km, at radius R from station
      
      h    = hh;
      for kk=1:nstat
        disp(sprintf('Station %4d: %8.2f  %8.2f',kk,stla(kk),stlo(kk)));
        % r2 = squared distance from station
        r2 = sqdist3(xx,yy,zz,xs(kk),ys(kk),zs(kk));
        tmp = hval0 + sqrt(r2)*(hvalR-hval0)/R;
        ii  = find(tmp<h); % min so far for requested length?
        h(ii) = tmp(ii);  % aber mit minus und find implementieren
      end
      % Copy to final matrix hh
      hh = h;
    end % if(linear_weighting)
    
    
    % Gaussian weighting was used for  grids USA2 to USA4. Was replaced
    % by linear weighting because transition from station-densified to
    % spherically symmetric area tended to be very abrupt.
    
    gauss_weighting = 0;
    if(gauss_weighting)  % good, use again later: Gaussian weighting functions 
      h    = zeros(size(hh));
      % Length scaling/resolving length (set above now)
      % (default: resolving length of upper mantle)
      %% sig = 200;    % width (sigma) of gausswin in km
      for kk=1:nstat
        disp(sprintf('Station %4d: %8.2f  %8.2f',kk,stla(kk),stlo(kk)));
        % r2 = squared distance from station
        r2 = sqdist3(xx,yy,zz,xs(kk),ys(kk),zs(kk));
        h = h + exp(-r2/(2*sig*sig)); % gaussian [0,1], width sig
      end
      % Threshold: if resolving length smaller than hmin then set to hmin.
      % If larger than local hh then set to hh
      % threshold
      fmx = 3; % maximum allowable densification factor
      h(h>fmx) = fmx;  
      h(h<1)   = 1;    % don't allow less dense grid

      % Add 3-D contribution of h to hh
      hh = hh./h;
    end % if(gauss_weighting)


    % End: resolution width function hh 

  %-----------------------
  case 'kerneldensified'
  %-----------------------
 
    disp('CAUTION for case = kerneldensified');
    disp('CAUTION: This option used to work and should work');
    disp('CAUTION: but I havent checked again. Also it seems');
    disp('CAUTION: to (over-)write to hardcoded file names -- FIX THAT!');
    error('Bailing out to be safe. Double-check and update code first.')
    
    
    % read in Guust's weighted kernel column sum
    meshdir = '/raid/AmplitudeProjects/testdistmesh/saved_meshes/USAmatrixdensified_0423';
    fres = 'reslength.new';

    ccdir = pwd; cd(meshdir);
    
    % either this...
    ffacet = 'facets.all';
    fparms = [fres '.parms'];
    fh     = [fres '.xyz'];
    fid    = fopen(fparms,'w');      
    fprintf(fid,'%s\n%s\n%d\n%d\n%s\n%d',...
	    fres,ffacet,100,1,fh,50);
    fclose(fid);
    disp(sprintf('\nInput for sol2xyz fortran program:'));
    disp(textread(fparms,'%s'));
    % Interpolate from tetrahedral to cubic grid
    fbin = '/data/jingjie/sigloch/INVERSION/Programs/buildgrid/sol2xyz';
    ucall = sprintf('cp %s %s.bkp; %s < %s',fh,fh,fbin,fparms); 
    disp(sprintf('\nCalling FORTRAN: \n%s\n',ucall));
    [flag,msg]= system(ucall)
    if(flag)
      error('ERROR: system call (FORTRAN) encountered error!');
    end
    hhh  = load(fh);  % read in again
    hdim = round(length(hhh)^(1/3)); % hhh = ndim x ndim x ndim
    hhh  = reshape(hhh,[hdim hdim hdim]);
    
    for k=1:hdim
      myx = [-6400:100:6400];
      tmp = squeeze(hhh(:,:,k));
      tmp = max(4*tmp,60);
      imagesc(tmp);
      title(k);
      colorbar vert
      pause
    end  
    
    if(0)
    % or that...
    fnam = fullfile(meshdir,fres);	    
    [X,Y,Z,myH] = textread(fnam,'%f%f%f%f','headerlines',2);
    hhh = griddata3(X,Y,Z,V,xx,yy,zz);  %% !! takes hours!
    save hhh.mat hhh xx yy zz X Y Z V 
    
    load hhh.mat
    end % if(0)
    
  %-----------------------
  otherwise
  %-----------------------
    error(sprintf('MODE NOT IMPLEMENTED: %s',hmatrix_mode))
  
end %  switch hmatrix_mode 
    

%--------------------------------------------
% Compute mesh points (distmesh algo) 
%--------------------------------------------


% h0: distance of initial random point distribution
% set equal to smallest resolving length
h0 = min(min(min(hh)));
	   
% MAKE MESH
disp(sprintf('Mesh gets saved every 5 iters. Okay to interrupt.'))
perm_vec = [2, 1, 3];
xx_p = permute(xx, perm_vec);
yy_p = permute(yy, perm_vec);
zz_p = permute(zz, perm_vec);

[p,t]=distmeshnd_hvar(d,@hmatrix3d,h0,...
	[-R0*ones(1,dim);R0*ones(1,dim)],...
        pfix,xx_p,yy_p,zz_p,hh,hh);
 
% Note: it may take extremely long or forever
% for distmeshnd_hvar to 
% reach its break condition. Instead, it can be interrupted
% whenever the mesh looks good enough (see plot that
% gives output every 5 iters). Every 5th iter also gets saved:
% save dmeshtmp p t

%-------------------------------------- 

% POST-PROCESSING


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
fh = 'hmatrix3d';  
HH0=feval(fh,p(:,:),xx_p,yy_p,zz_p,hh,hh);

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

 
 
 
 
