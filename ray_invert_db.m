function ray_invert_db(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_invert_db
%
% This function is a front-end to the raytrace3d invert routine using arrivals
% from a css3.0 database. Note: this routine currently only extracts and 
% inverts for P arrivals (although it can be easily modified to include S 
% waves as well). This function requires Antelope to be installed on the
% system.
%
% This function requires the presence of two other functions:
%   ray_make_traveltime.m
%   ray_latlon2xyz_flat.m (for flat-earth projection)
%       OR
%   ray_latlon2xyz.m (for spherical-earth projection)
%
% Usage:
% ray_invert_db(database,model_file,model_pars_file,damping,[subset],[locations])
%
% Required Inputs:
%   database:        The path to a css3.0 database containing the following
%                    tables: site, affiliation, arrival, assoc, and origin
%
%   model_file:      The starting 3D model file used in any raytrace3d routine
%
%   model_pars_file: The path to the input_mps.txt file required by
%                    raytrace3d invert (see raytrace3d help for file
%                    description)
%
%   damping:         The value of the damping parameter to be used in the
%                    inversion
%
% Optional Input:
%   subset:          A valid string for subsetting a Datascope database
%
%   locations:       The path to the output file containing updated
%                    earthquake hypocenters generated by the ray_locate 
%                    function
%
% Output (to filesystem):
%   inversion_model.model:      The model result from the inversion
%
%   inversion_model.model.hits: The number of ray hits at each model node
%
%
% Author:
% Matt Gardine
% February 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks that Antelope is installed on the system.
if (exist('dbopen') ~= 3)
    error('Error: Antelope must be installed on the system to use this function')
end

% Checks for the existence of the ray_make_traveltimes function 
if (exist('ray_make_traveltime') ~= 2)
    error('Error: This function is dependent on ray_make_traveltime. Please add this function into the path')
end

% Checks for the existence of the ray_defaults file
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

% Checks for the existence of ray_latlon2xyz or ray_latlon2xyz_flat
if strcmp(projection,'flat')
    if (exist('ray_latlon2xyz_flat') ~= 2)
        error('Error: This function is dependent on ray_latlon2xyz_flat.  Please add this function into the path')
    end
elseif strcmp(projection,'spherical')
    if (exist('ray_latlon2xyz') ~= 2)
        error('Error: This function is dependent on ray_latlon2xyz.  Please add this function into the path')
    end
end

switch nargin
    case 4
        database=varargin{1};
        model_file=varargin{2};
        model_pars_file=varargin{3};
        damping=varargin{4};
        
        db = dbopen(database,'r');
        db1a = dblookup(db,'','site','','');
        db1b = dblookup(db,'','affiliation','','');
        db1c = dblookup(db,'','arrival','','');
        db1d = dblookup(db,'','assoc','','');
        db1e = dblookup(db,'','origin','','');

        db = dbjoin(db1a,db1b);
        db = dbjoin(db,db1c);
        db = dbjoin(db,db1d);
        db = dbjoin(db,db1e);

        db = dbsubset(db,'phase=~/P/');
        
        ray_make_traveltime(db);
        
    case 5
        database=varargin{1};
        model_file=varargin{2};
        model_pars_file=varargin{3};
        damping=varargin{4};
        
        db = dbopen(database,'r');
        db1a = dblookup(db,'','site','','');
        db1b = dblookup(db,'','affiliation','','');
        db1c = dblookup(db,'','arrival','','');
        db1d = dblookup(db,'','assoc','','');
        db1e = dblookup(db,'','origin','','');

        db = dbjoin(db1a,db1b);
        db = dbjoin(db,db1c);
        db = dbjoin(db,db1d);
        db = dbjoin(db,db1e);

        db = dbsubset(db,'phase=~/P/');
        
        if exist(varargin{5},'file')~=2
            subset=varargin{5};
            db = dbsubset(db,subset);
            ray_make_traveltime(db);
        else
            locations=varargin{5};
            ray_make_traveltime(db,locations);
        end
        
    case 6
        database=varargin{1};
        model_file=varargin{2};
        model_pars_file=varargin{3};
        damping=varargin{4};
        
        if exist(varargin{5},'file')~=2
            subset=varargin{5};
            locations=varargin{6};
        else
            locations=varargin{5};
            subset=varargin{6};
        end
        
        db = dbopen(database,'r');
        db1a = dblookup(db,'','site','','');
        db1b = dblookup(db,'','affiliation','','');
        db1c = dblookup(db,'','arrival','','');
        db1d = dblookup(db,'','assoc','','');
        db1e = dblookup(db,'','origin','','');

        db = dbjoin(db1a,db1b);
        db = dbjoin(db,db1c);
        db = dbjoin(db,db1d);
        db = dbjoin(db,db1e);

        db = dbsubset(db,'phase=~/P/');
        db = dbsubset(db,subset);
        
        ray_make_traveltime(db,locations);
end

runstring = ['raytrace3d invert ' model_file ' ' model_pars_file ' ./traveltimes.tsv 5 ' num2str(damping)...
    ' 0 180 100 0 360 181 ./inversion_model.model >& invert_output.log'];
system(runstring);
