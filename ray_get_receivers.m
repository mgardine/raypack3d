function ray_get_receivers(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_get_receivers
%
% This function is internal to the ray_shoot_predicted function. It generates
% the receiver.rec file as required by raytrace3d source_to_receiver.
%
% This function requires the presence of one other function:
%   ray_latlon2xyz_flat.m (for flat-earth projection)
%       OR
%   ray_latlon2xyz.m (for spherical-earth projection)
%
% Author: 
% Matt Gardine
% August 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('ray_defaults','file')==2
    [ref_lat,ref_lon,projection]=ray_defaults();
else
    ref_lat=17.01;
    ref_lon=-105.99;
    projection='flat';
end

switch nargin
    case 5
        database = varargin{1};
        orid = varargin{2};
        outfile = varargin{3};
        ref_lat = varargin{4};
        ref_lon = varargin{5};

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

        subset = ['orid=~/' num2str(orid) '/ && phase=~/P/'];
        
        db = dbsubset(db,subset);

        [lat,lon,elev,orig_time,arr_time] = dbgetv(db,'lat','lon','elev','origin.time','time');
        dbclose(db);
        
    case 6
        database = varargin{1};
        orid = varargin{2};
        outfile = varargin{3};
        ref_lat = varargin{4};
        ref_lon = varargin{5};
        subset = varargin{6};

        db = dbopen(database,'r');
        
        db1a = dblookup(db,'','origin','','');
        db1b = dblookup(db,'','assoc','','');
        db1c = dblookup(db,'','arrival','','');
        db1d = dblookup(db,'','site','','');
        db1e = dblookup(db,'','affiliation','','');

        db = dbjoin(db1a,db1b);
        db = dbjoin(db,db1c);
        db = dbjoin(db,db1d);
        db = dbjoin(db,db1e);

        subset = ['orid=~/' num2str(orid) '/ && phase=~/P/ && ' subset ];
        
        db = dbsubset(db,subset);

        [lat,lon,elev,orig_time,arr_time] = dbgetv(db,'site.lat','site.lon','elev','origin.time','arrival.time');
        dbclose(db);
        
    otherwise
        disp('Invalid options for get_receivers');
end

if strcmp(projection,'flat')
    fid = fopen(outfile,'wt');
    fprintf(fid,'%s\n','x     y     z     T0     T     ratio     order');
    for i=1:length(lat)
        [x,y,z] = ray_latlon2xyz_flat(lat(i),lon(i),elev(i),ref_lat,ref_lon);
        fprintf(fid,'%f %f %f %f %f 1 0\n',x,y,z,orig_time(i),arr_time(i));
    end
    
elseif strcmp(projection,'spherical')
fid = fopen(outfile,'wt');
    fprintf(fid,'%s\n','x     y     z     T0     T     ratio     order');
    for i=1:length(lat)
        [x,y,z] = ray_latlon2xyz(lat(i),lon(i),elev(i),ref_lat,ref_lon);
        fprintf(fid,'%f %f %f %f %f 1 0\n',x,y,z,orig_time(i),arr_time(i));
    end
    
else
    disp('Error: Invalid projection type');
    return
end

fclose(fid);