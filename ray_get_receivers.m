function ray_get_receivers(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_get_receivers
%
% This function is internal to the ray_shoot_predicted function. It generates
% the receiver.rec file as required by raytrace3d source_to_receiver.
%
% This function requires the presence of one other function:
%   ray_latlon2xyz.m
%
% Author: 
% Matt Gardine
% August 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 7
        database = varargin{1};
        orid = varargin{2};
        outfile = varargin{3};
        ref_lat = varargin{4};
        ref_lon = varargin{5};
        ref_alt = varargin{6};
        projection = varargin{7};

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
        
    case 8
        database = varargin{1};
        orid = varargin{2};
        outfile = varargin{3};
        ref_lat = varargin{4};
        ref_lon = varargin{5};
        ref_alt = varargin{6};
        projection = varargin{7};
        subset = varargin{8};

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

   

fid = fopen(outfile,'wt');
fprintf(fid,'%s\n','x     y     z     T0     T     ratio     order');
for i=1:length(lat)
    [x,y,z] = ray_latlon2xyz(lat(i),lon(i),elev(i),ref_lat,ref_lon,ref_alt,projection);
    fprintf(fid,'%f %f %f %f %f 1 0\n',x,y,z,orig_time(i),arr_time(i));
end

fclose(fid);