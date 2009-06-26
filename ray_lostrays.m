function lostray=ray_lostrays(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_lostrays
%
% This function tracks the origin-station pairs in the original database 
% that are unable to be matched by the raytrace routine.
%
% This function requires the presence of one other function:
%   ray_latlon2xyz_flat.m (for flat-earth projection)
%       OR
%   ray_latlon2xyz.m (for spherical-earth projection)
%
% Usage: lostray=ray_lostrays(r,database)
%
% Required Inputs:
%   r:              A raytrace3d structure with the hypocenter and station
%                   arrays filled.
%
%   database:       Full path to the css3.0 database.  This database must
%                   have working origin, site, affiliation, assoc, and 
%                   arrival tables.
%
%
% Output:
%   lostray:        An array of origin-station pairs in the database that 
%                   are unable to be matched by the raytrace routine.
%                   Format: 
%                   [orig_x orig_y orig_z station_x station_y orig_time]
%
%
% Author:
% Matt Gardine
% April 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('ray_defaults','file')==2
    [ref_lat,ref_lon,projection]=ray_defaults();
    disp('ray_defaults file found.')
    disp(['ref_lat = ' num2str(ref_lat)])
    disp(['ref_lon = ' num2str(ref_lon)])
    disp(['projection = ' projection])
else
    ref_lat=17.01;
    ref_lon=-105.99;
    projection='flat';
    disp('ray_defaults file NOT found. Values used:')
    disp(['ref_lat = ' num2str(ref_lat)])
    disp(['ref_lon = ' num2str(ref_lon)])
    disp(['projection = ' projection])
end


switch nargin
    case 2
        r=varargin{1};
        database=varargin{2};
        ref_lat=17.01;
        ref_lon=-105.99;
    case 4
        r=varargin{1};
        database=varargin{2};
        ref_lat=varargin{3};
        ref_lon=varargin{4};
    otherwise
        disp('Invalid number in inputs. Usage: lostray=ray_lostrays(r.database,[ref_lat],[ref_lon]')
end

db = dbopen(database,'r');
db1a = dblookup(db,'','site','','');
db1b = dblookup(db,'','affiliation','','');
db1c = dblookup(db,'','origin','','');
db1d = dblookup(db,'','assoc','','');
db1e = dblookup(db,'','arrival','','');

db = dbjoin(db1a,db1b);
db = dbjoin(db,db1e);
db = dbjoin(db,db1d);
db = dbjoin(db,db1c);

dbsubset(db,'phase=~/P/');
db = dbsort(db,'orid');

[orig_lat,orig_lon,depth,site_lat,site_lon,orig_time] = dbgetv(db,'origin.lat','origin.lon','depth','site.lat','site.lon','origin.time');
dbclose(db);

if strcmp(projection,'flat')
    for i=1:length(orig_lat)
        [orig_x,orig_y,junk]=ray_latlon2xyz_flat(orig_lat(i),orig_lon(i),depth(i),ref_lat,ref_lon);
        [site_x,site_y,junk]=ray_latlon2xyz_flat(site_lat(i),site_lon(i),0,ref_lat,ref_lon);
        dbase(i,:)=[orig_x orig_y depth site_x site_y orig_time(i)];
    end
    
elseif strcmp(projection,'spherical')
    for i=1:length(orig_lat)
        [orig_x,orig_y,junk]=ray_latlon2xyz(orig_lat(i),orig_lon(i),depth(i),ref_lat,ref_lon);
        [site_x,site_y,junk]=ray_latlon2xyz(site_lat(i),site_lon(i),0,ref_lat,ref_lon);
        dbase(i,:)=[orig_x orig_y depth site_x site_y orig_time(i)];
    end
    
else
    disp('Error: Invalid projection type');
    return
end

i=1;
k=1;
while i<length(r.hypocenter(:,1))
    a = find(r.hypocenter(:,1)<=(r.hypocenter(i,1)+0.0001) & r.hypocenter(:,1)>=(r.hypocenter(i,1)-0.0001) ...
        & r.hypocenter(:,2)<=(r.hypocenter(i,2)+0.0001) & r.hypocenter(:,2)>=(r.hypocenter(i,2)-0.0001));
    b = find(dbase(:,1)<=(r.hypocenter(i,1)+0.0001) & dbase(:,1)>=(r.hypocenter(i,1)-0.0001) ...
        & dbase(:,2)<=(r.hypocenter(i,2)+0.0001) & dbase(:,2)>=(r.hypocenter(i,2)-0.0001));
    pred_stations=[r.station(a,1) r.station(a,2)];
    db_origins=[dbase(b,1) dbase(b,2) dbase(b,3) dbase(b,4) dbase(b,5) dbase(b,6)];
    for j=1:size(db_origins,1)
        c = pred_stations(:,1)<=(db_origins(j,4)+0.0001) & pred_stations(:,1)>=(db_origins(j,4)-0.0001) ...
            & pred_stations(:,2)<=(db_origins(j,5)+0.0001) & pred_stations(:,2)>=(db_origins(j,5)-0.0001);
        if ~any(c)
            lostray(k,:)=[db_origins(j,1) db_origins(j,2) db_origins(j,3) db_origins(j,4) db_origins(j,5) db_origins(j,6)];
            k=k+1;
        end
    end
    i=i+length(a);
    clear a b pred_stations db_origins
end
