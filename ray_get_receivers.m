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

        subset = ['orid=~/' num2str(orid) '/ && phase=~/P/'];
        
        db = dbsubset(db,subset);

        [lat,lon,elev,orig_time,arr_time] = dbgetv(db,'site.lat','site.lon','site.elev','origin.time','arrival.time');
        dbclose(db);
        
    case 6
        database = varargin{1};
        orid = varargin{2};
        outfile = varargin{3};
        ref_lat = varargin{4};
        ref_lon = varargin{5};
        subset = varargin{6};

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

        subset = ['orid=~/' num2str(orid) '/ && phase=~/P/ && ' subset ];
        
        db = dbsubset(db,subset);

        [lat,lon,elev,orig_time,arr_time] = dbgetv(db,'site.lat','site.lon','site.elev','origin.time','arrival.time');
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