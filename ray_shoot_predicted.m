function ray_shoot_predicted(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_shoot_predicted
%
% This function shoots rays from origins to stations in a css3.0 database
% using the raytrace3d shoot_star and raytrace3d source_to_receiver
% functions. It generates four files as part of the raytrace3d suite:
% receiver_file.rec, star_file.star, pred_tt.tt, ray_file.rays.
%
% This function requires the presence of two other functions:
%   ray_get_receivers.m
%   ray_latlon2xyz_flat.m (for flat-earth projection)
%       OR
%   ray_latlon2xyz.m (for spherical-earth projection)
% 
% Usage: ray_shoot_predicted(database,model_file,[station])
%
% Inputs:
%   database:       Full path to the css3.0 database.  This database must have
%                   working origin, site, affiliation, assoc, and arrival tables.
%
%   model_file:     A 3-D velocity model file in the proper format required
%                   by raytrace3d shoot_star (see raytrace3d shoot_star for
%                   format)
%
%
% Optional Inputs:
%   station:        A station name in the css3.0 site table to subset on. 
%
%
% Output:
%   origin_data:    A 4 column array, with origin latitides, longitudes, times,
%                   and magnitudes.
%
%
% Author:
% Matt Gardine
% August 2008
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
        database = varargin{1};
        model_file = varargin{2};
        
        db = dbopen(database,'r');
        db = dblookup(db,'','origin','','');

        [lat,lon,depth,orids] = dbgetv(db,'lat','lon','depth','orid');
        dbclose(db);
        
        ttout = fopen('./pred_tt.tt','wt');
        rayout = fopen('./ray_file.rays','wt');

        fprintf(ttout,'%s\n','id    x0    y0    z0    x1    y1    z1    T0    Ta    TT    T    dT    ratio   obs_order    pre_order    amplitude   ntts');
        fprintf(rayout,'%s\n','id    x    y    z    tetra');

        for i=1:length(orids)
            disp(['orid is: ' num2str(orids(i))]);
            ray_get_receivers(database,orids(i),'./receiver_file.rec',ref_lat,ref_lon);
            
            if strcmp(projection,'flat')
                [x,y,z]=ray_latlon2xyz_flat(lat(i),lon(i),-1*depth(i),ref_lat,ref_lon);
                
            elseif strcmp(projection,'spherical')
                [x,y,z]=ray_latlon2xyz(lat(i),lon(i),-1*depth(i),ref_lat,ref_lon);
                
            else
                disp('Error: Invalid projection type');
            end
            
            runstring = ['raytrace3d shoot_star ' model_file ' ' num2str(x) ' '...
                num2str(y) ' ' num2str(z) ' 0 180 181 0 360 361 0 > ./star_file.star'];
            system(runstring);
            runstring2 = ['raytrace3d source_to_receivers ' model_file ' ' num2str(x) ' '...
                num2str(y) ' ' num2str(z) ' ./star_file.star 0 ./receiver_file.rec ' ...
                './pred_tt_temp.tt ./ray_file_temp.rays'];
            system(runstring2);
            [idr Xr Yr Zr tetrar ]=textread('./ray_file_temp.rays','%d %f %f %f %d','headerlines',1);
            [id x0 y0 z0 x1 y1 z1 T0 Ta TT T dT ratio obs_order pre_order amplitude ntts]=textread('./pred_tt_temp.tt','%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d','headerlines',1);
            system('rm -f ./ray_file_temp.rays ./pred_tt_temp.tt ./star_file.star ./receiver_file.rec');
            for j=1:length(idr)
                fprintf(rayout,'%d %f %f %f %d\n',idr(j),Xr(j),Yr(j),Zr(j),tetrar(j));
            end
            for j=1:length(id)
                fprintf(ttout,'%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d\n',id(j),x0(j),y0(j),z0(j),x1(j),y1(j),z1(j),T0(j),Ta(j),TT(j),T(j),dT(j),ratio(j),obs_order(j),pre_order(j),amplitude(j),ntts(j));
            end
        end
        fclose(rayout);
        fclose(ttout);
        
    case 3
        database = varargin{1};
        model_file = varargin{2};
        station = varargin{3};
        
        db = dbopen(database,'r');
        
        db2a = dblookup(db,'','origin','','');
        db2b = dblookup(db,'','assoc','','');
        db2c = dblookup(db,'','arrival','','');
        db2d = dblookup(db,'','site','','');
        db2e = dblookup(db,'','affiliation','','');

        db2 = dbjoin(db2a,db2b);
        db2 = dbjoin(db2,db2c);
        db2 = dbjoin(db2,db2d);
        db2 = dbjoin(db2,db2e);
        
        subset = ['phase=~/P/&&sta=~/' station '/'];
        db2 = dbsubset(db2,subset);
        
        [lat,lon,depth,orids] = dbgetv(db2,'lat','lon','depth','orid');
        dbclose(db);
        
        ttout = fopen('./pred_tt.tt','wt');
        rayout = fopen('./ray_file.rays','wt');

        fprintf(ttout,'%s\n','id    x0    y0    z0    x1    y1    z1    T0    Ta    TT    T    dT    ratio   obs_order    pre_order    amplitude   ntts');
        fprintf(rayout,'%s\n','id    x    y    z    tetra');

        for i=1:length(orids)
            disp(['orid is: ' num2str(orids(i))]);
            ray_get_receivers(database,orids(i),'./receiver_file.rec',ref_lat,ref_lon,station);
            
            if strcmp(projection,'flat')
                [x,y,z]=ray_latlon2xyz_flat(lat(i),lon(i),-1*depth(i),ref_lat,ref_lon);
                
            elseif strcmp(projection,'spherical')
                [x,y,z]=ray_latlon2xyz(lat(i),lon(i),-1*depth(i),ref_lat,ref_lon);
                
            else
                disp('Error: Invalid projection type');
            end
            
            runstring = ['raytrace3d shoot_star ' model_file ' ' num2str(x) ' '...
                num2str(y) ' ' num2str(z) ' 0 180 181 0 360 361 0 > ./star_file.star'];
            system(runstring);
            runstring2 = ['raytrace3d source_to_receivers ' model_file ' ' num2str(x) ' '...
                num2str(y) ' ' num2str(z) ' ./star_file.star 0 ./receiver_file.rec ' ...
                './pred_tt_temp.tt ./ray_file_temp.rays'];
            system(runstring2);
            [idr Xr Yr Zr tetrar ]=textread('./ray_file_temp.rays','%d %f %f %f %d','headerlines',1);
            [id x0 y0 z0 x1 y1 z1 T0 Ta TT T dT ratio obs_order pre_order amplitude ntts]=textread('./pred_tt_temp.tt','%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d','headerlines',1);
            system('rm -f ./ray_file_temp.rays ./pred_tt_temp.tt ./raytrace3d/star_file.star ./receiver_file.rec');
            for j=1:length(idr)
                fprintf(rayout,'%d %f %f %f %d\n',idr(j),Xr(j),Yr(j),Zr(j),tetrar(j));
            end
            for j=1:length(id)
                fprintf(ttout,'%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d\n',id(j),x0(j),y0(j),z0(j),x1(j),y1(j),z1(j),T0(j),Ta(j),TT(j),T(j),dT(j),ratio(j),obs_order(j),pre_order(j),amplitude(j),ntts(j));
            end
        end
        fclose(rayout);
        fclose(ttout);
        
    otherwise
        disp('Invalid inputs')
        disp('Usage: ray_shoot_predicted(database,model_file,[station]')
end
