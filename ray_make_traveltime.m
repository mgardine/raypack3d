function ray_make_traveltime(db)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_make_traveltime
%
% This function is internal to the ray_invert_db and ray_locate functions.
% It generates the traveltimes.tsv file required by raytrace3d invert.
%
% This function requires the presence of one other function:
%   ray_latlon2xyz_flat.m (for flat-earth projection)
%       OR
%   ray_latlon2xyz.m (for spherical-earth projection)
%
% Input:
%   db:         A database pointer containing, at a minimum, the site, 
%               arrival, assoc, and origin tables.
%
% Author: 
% Matt Gardine
% June 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('ray_defaults','file')==2
    [ref_lat,ref_lon,projection]=ray_defaults();
else
    ref_lat=17.01;
    ref_lon=-105.99;
    projection='flat';
end

db = dbsort(db,'orid');

[orig_lat,orig_lon,orig_depth,sta_lat,sta_lon,sta_elev,orig_time,arr_time]=...
    dbgetv(db,'origin.lat','origin.lon','depth','lat','lon','elev','origin.time','time');

fid = fopen('./traveltimes.tsv','wt');

orig_x=zeros(length(orig_lat),1);
orig_y=zeros(length(orig_lat),1);
orig_z=zeros(length(orig_lat),1);
sta_x=zeros(length(orig_lat),1);
sta_y=zeros(length(orig_lat),1);
sta_z=zeros(length(orig_lat),1);
time=zeros(length(orig_lat),1);

if strcmp(projection,'flat')
    for i=1:length(orig_lat)
        [orig_x(i),orig_y(i),orig_z(i)]=ray_latlon2xyz_flat(orig_lat(i),orig_lon(i),-1*orig_depth(i),ref_lat,ref_lon);
        [sta_x(i),sta_y(i),sta_z(i)]=ray_latlon2xyz_flat(sta_lat(i),sta_lon(i),sta_elev(i),ref_lat,ref_lon);
        time(i)=arr_time(i)-orig_time(i);
        fprintf(fid,'%f %f %f %f %f %f 0 0 %f 1 0\n',orig_x(i),orig_y(i),orig_z(i),sta_x(i),sta_y(i),sta_z(i),time(i));
    end
    
elseif strcmp(projection,'spherical')
    for i=1:length(orig_lat)
        [orig_x(i),orig_y(i),orig_z(i)]=ray_latlon2xyz(orig_lat(i),orig_lon(i),-1*orig_depth(i),ref_lat,ref_lon);
        [sta_x(i),sta_y(i),sta_z(i)]=ray_latlon2xyz(sta_lat(i),sta_lon(i),sta_elev(i),ref_lat,ref_lon);
        time(i)=arr_time(i)-orig_time(i);
        fprintf(fid,'%f %f %f %f %f %f 0 0 %f 1 0\n',orig_x(i),orig_y(i),orig_z(i),sta_x(i),sta_y(i),sta_z(i),time(i));
    end
    
else
    disp('Error: Invalid projection type');
end

fclose(fid);