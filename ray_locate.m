function ray_locate(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_locate
%
% This function locates origins using arrival data found in a css3.0
% database and a raytrace3d model file.
%
% This function requires the presence of one other functions:
%   ray_latlon2xyz.m
%
% Usage:
%   ray_locate(database,model_file,tolerance,damping,[subset])
%
% Required Inputs:
%   database:         The path to a css3.0 database containing the
%                           following: tables: site, affiliation, arrival, 
%                           assoc, origin, and event
%
%   model_file:             The 3D model file used in any raytrace3d routine
%
%   tolerance:         The value, in seconds, for which the iterating
%                           will stop when the rms traveltime error is less 
%                           than.
%
%   damping:         A fraction of the maximum GTG value.
%
% Optional Input:
%   subset:                 A valid string for subsetting a Datascope database
%
% Output (to filesystem):
%   location_output.tsv:    The new origin data
%
% Author:
% Matt Gardine
% June 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks that Antelope is installed on the system.
if (exist('dbopen') ~= 3)
    error('Error: Antelope must be installed on the system to use this function')
end

% Checks for the existence of the ray_defaults file
if exist('ray_defaults','file')==2
    [ref_lat,ref_lon,ref_alt,projection]=ray_defaults();
    disp('ray_defaults file found.')
    disp(['ref_lat = ' num2str(ref_lat)])
    disp(['ref_lon = ' num2str(ref_lon)])
    disp(['ref_alt = ' num2str(ref_alt)])
    disp(['projection = ' projection])
else
    ref_lat=17.01;
    ref_lon=-105.99;
    ref_alt=0;
    projection='flat';
    disp('ray_defaults file NOT found. Values used:')
    disp(['ref_lat = ' num2str(ref_lat)])
    disp(['ref_lon = ' num2str(ref_lon)])
    disp(['ref_alt = ' num2str(ref_alt)])
    disp(['projection = ' projection])
end

% Checks for the existence of ray_latlon2xyz
if (exist('ray_latlon2xyz') ~= 2)
    error('Error: This function is dependent on ray_latlon2xyz.  Please add this function into the path')
end

switch nargin
    case 4
        database=varargin{1};
        model_file=varargin{2};
        tolerance=varargin{3};
        damping=varargin{4};
        
        db = dbopen(database,'r');
        db_site = dblookup(db,'','site','','');
        db_affil = dblookup(db,'','affiliation','','');
        db_arr = dblookup(db,'','arrival','','');
        db_assoc = dblookup(db,'','assoc','','');
        db_orig = dblookup(db,'','origin','','');
        db_event = dblookup(db,'','event','','');
        
        db = dbjoin(db_orig,db_event);
        
        db = dbsubset(db,'orid==prefor');
        
        db = dbjoin(db,db_assoc);
        db = dbjoin(db,db_arr);
        db = dbjoin(db,db_affil);
        db = dbjoin(db,db_site);
           
    case 5
        database=varargin{1};
        model_file=varargin{2};
        tolerance=varargin{3};
        damping=varargin{4};
        subset=varargin{5};
        
        db = dbopen(database,'r');
        db_site = dblookup(db,'','site','','');
        db_affil = dblookup(db,'','affiliation','','');
        db_arr = dblookup(db,'','arrival','','');
        db_assoc = dblookup(db,'','assoc','','');
        db_orig = dblookup(db,'','origin','','');
        db_event = dblookup(db,'','event','','');
        
        db = dbjoin(db_orig,db_event);
        
        db = dbsubset(db,'orid==prefor');
        
        db = dbjoin(db,db_assoc);
        db = dbjoin(db,db_arr);
        db = dbjoin(db,db_affil);
        db = dbjoin(db,db_site);

        db = dbsubset(db,subset);       
        
    otherwise
        help ray_locate
        return;
end

orids=dbgetv(db,'orid');
orids=unique(orids);

fout=fopen('./location_output.tsv','wt');
fprintf(fout,'%s\n','orid xs ys zs t0 sigma-x sigma-y sigma-z sigma-t n-data iter rms');

for i=1:length(orids)
    db_temp=dbsubset(db,['orid=~/' num2str(orids(i)) '/']);
    
    [orig_lat,orig_lon,orig_depth,orig_time]=dbgetv(db_temp,'origin.lat','origin.lon','origin.depth','origin.time');
        
    [xs,ys,zs]=ray_latlon2xyz(orig_lat(1),orig_lon(1),-1*orig_depth(1),ref_lat,ref_lon,ref_alt,projection);
    
    make_ttfile(db_temp)
            
    runstring = ['raytrace3d locate ' model_file ' ' num2str(xs) ' ' num2str(ys) ' ' num2str(zs) ' ' ...
    ' 0 180 100 0 360 181 ./temp_traveltimes.tsv ' num2str(tolerance) ' ./temp_location.tsv ' num2str(damping) ...
    ' 1 0 5' ];
    system(runstring);
    
    [xs ys zs t0 sigx sigy sigz sigt n iter rms]=textread('./temp_location.tsv','%f %f %f %f %f %f %f %f %d %d %f','headerlines',1);
    fprintf(fout,'%d %f %f %f %f %f %f %f %f %d %d %f\n',orids(i),xs,ys,zs,orig_time(1)+t0,sigx,sigy,sigz,sigt,n,iter,rms);
    
    system('rm -f temp_traveltimes.tsv temp_location.tsv');
    
end

fclose(fout);
dbclose(db);

function make_ttfile(db)

if exist('ray_defaults','file')==2
    [ref_lat,ref_lon,projection]=ray_defaults();
else
    ref_lat=17.01;
    ref_lon=-105.99;
    projection='flat';
end

[sta_lat,sta_lon,sta_elev,orig_time,arr_time,phase]=dbgetv(db,'site.lat','site.lon','site.elev','origin.time','arrival.time','phase');

fid = fopen('./temp_traveltimes.tsv','wt');
fprintf(fid,'%s\n','x y z T0 T ratio order');

for i=1:length(sta_lat)
    [sta_x(i),sta_y(i),sta_z(i)]=ray_latlon2xyz(sta_lat(i),sta_lon(i),sta_elev(i),ref_lat,ref_lon,ref_alt,projection);
    time(i)=arr_time(i)-orig_time(i);
    if strcmp(phase(i),'P')
        fprintf(fid,'%f %f %f 0 %f 1 0\n',sta_x(i),sta_y(i),sta_z(i),time(i));
    elseif strcmp(phase(i),'S')
        fprintf(fid,'%f %f %f 0 %f 1.76 0\n',sta_x(i),sta_y(i),sta_z(i),time(i));
    end
end


fclose(fid);