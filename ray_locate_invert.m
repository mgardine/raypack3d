function ray_locate_invert(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_locate_invert
%
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

% Checks for the existence of the ray_make_traveltimes function 
if (exist('ray_make_traveltime') ~= 2)
    error('Error: This function is dependent on ray_make_traveltime. Please add this function into the path')
end

switch nargin
    case 5
        database=varargin{1};
        model_file=varargin{2};
        model_pars_file=varargin{3};
        damping=varargin{4};
        rel_weight=varargin{5};
        
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
        
        event_from_db=1;
              
    case 6
        database=varargin{1};
        model_file=varargin{2};
        model_pars_file=varargin{3};
        damping=varargin{4};
        rel_weight=varargin{5};
        
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
        
        if exist(varargin{6},'file')~=2
            subset=varargin{6};
            disp(['Subsetting database using expression "' subset '"']);
            db = dbsubset(db,subset);
            event_from_db=1;
        else
            locations=varargin{6};
            disp(['Using events file ' locations]);
            event_from_db=0;
        end
        
    case 7
        database=varargin{1};
        model_file=varargin{2};
        model_pars_file=varargin{3};
        damping=varargin{4};
        rel_weight=varargin{5};
        
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
        
        locations=varargin{6};
        subset=varargin{7};
        
        disp(['Subsetting database using expression "' subset '"']);
        disp(['Using events file ' locations]);
        
        event_from_db=0;
        db = dbsubset(db,subset);
end

db=dbsort(db,'sta');
[tmp_sta,tmp_lat,tmp_lon,tmp_elev]=dbgetv(db,'sta','site.lat','site.lon','site.elev');
[sta,m,n]=unique(tmp_sta);
sta_lat=tmp_lat(m);
sta_lon=tmp_lon(m);
sta_elev=tmp_elev(m);

fid_station = fopen('./station.txt','wt');

for i=1:length(sta_lat)
   disp(['Working on station ' char(sta(i))])
   
   [sta_x,sta_y,sta_z]=ray_latlon2xyz(sta_lat(i),sta_lon(i),sta_elev(i),ref_lat,ref_lon,ref_alt,projection);
   
%   sta_z=sta_z+0.5;
   disp('Making tt_table')
   system(['(/home/menke/raytrace3d/RAYTRACE3D/raytrace3d tt_table ' model_file ' ' num2str(sta_x) ...
      ' ' num2str(sta_y) ' ' num2str(sta_z) ' 0 180 100 0 360 181 0.5 120 > tt_table_' char(sta(i)) ') >& tt_table_' char(sta(i)) '.log']);
   
%    disp('Making tt_index')
%    system(['/home/menke/raytrace3d/RAYTRACE3D/raytrace3d index_traveltimes ' model_file ' tt_table tt_index']);
   
   disp('Making jhd_table')   
   if event_from_db
      db_temp=dbsubset(db,['sta=~/' char(sta(i)) '/ && phase=~/P/']);
      [tmp_orid,tmp_lat,tmp_lon,tmp_depth]=dbgetv(db_temp,'origin.orid','origin.lat','origin.lon','origin.depth');
   
      [orid,m,n]=unique(tmp_orid);
      lat=tmp_lat(m);
      lon=tmp_lon(m);
      depth=tmp_depth(m); 
        
      fid = fopen('./events_temp.txt','wt');
      for j=1:length(orid)
          [orig_x,orig_y,orig_z]=ray_latlon2xyz(lat(j),lon(j),-1*depth(j),ref_lat,ref_lon,ref_alt,projection);
          fprintf(fid,'%d %.2f %.2f %.2f 0 120\n',orid(j),orig_x,orig_y,orig_z);
      end
      fclose(fid);
   else
      [orid orig_x orig_y orig_z time]=textread(locations,'%d %f %f %f %f');
      fid = fopen('./events_temp.txt','wt');
      for j=1:length(orid)
          fprintf(fid,'%d %.2f %.2f %.2f %f 120\n',orid(j),orig_x(j),orig_y(j),orig_z(j),time(j));
      end
      fclose(fid);
   end
   
   system(['(/home/menke/raytrace3d/RAYTRACE3D/raytrace3d jhd_table ' model_file ' tt_table_' char(sta(i)) ' e < events_temp.txt > jhd_table_' char(sta(i)) '.txt) >& jhd_table_' char(sta(i)) '.log']);
   
   system('rm -rf events_temp.txt');
   
   fprintf(fid_station,'%d %.2f %.2f %.2f %s\n',i,sta_x,sta_y,sta_z,['jhd_table_' char(sta(i)) '.txt']);
   
end

fclose(fid_station);

if event_from_db
    [tmp_orid,tmp_lat,tmp_lon,tmp_depth,tmp_time]=dbgetv(db,'orid','origin.lat','origin.lon','origin.depth','origin.time');
    [orid,m,n]=unique(tmp_orid);
    lat=tmp_lat(m);
    lon=tmp_lon(m);
    depth=tmp_depth(m);
    time=tmp_time(m);
    
    fid = fopen('./events.txt','wt');
    for j=1:length(orid)
        [orig_x,orig_y,orig_z]=ray_latlon2xyz(lat(j),lon(j),-1*depth(j),ref_lat,ref_lon,ref_alt,projection);
        fprintf(fid,'%d %.2f %.2f %.2f %.2f \n',orid(j),orig_x,orig_y,orig_z,0);
    end
    fclose(fid);
    locations='./events.txt';
end

% Old way of making phase.txt

% fid_phase = fopen('./phase.txt','wt');
% 
% db=dbsubset(db,'phase=~/P/');
% db=dbsort(db,'sta','orid');
% [sta_tmp,orid,orig_time,arr_time]=dbgetv(db,'sta','orid','origin.time','arrival.time');
% 
% j=1;
% for i=1:length(sta)
%     while (j<=length(orid) && strcmp(char(sta_tmp(j)),char(sta(i))))
%         fprintf(fid_phase,'1 %.2f 10 10\n',arr_time(j)-orig_time(j));
%         fprintf(fid_phase,'%d %d 0 1.0 1.0\n',i,orid(j));
%         j=j+1;
%     end
% end
% 
% fclose(fid_phase);

fid_phase = fopen('./phase.txt','wt');

db=dbsubset(db,'phase=~/P/');
db=dbsort(db,'sta','orid');

[ev_orid,ev_time]=textread(locations,'%d %*f %*f %*f %f');

for i=1:length(sta)
    [tmpid]=textread(['jhd_table_' char(sta(i)) '.txt'],'%d %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f');
    tmpid=unique(tmpid);
    sta_sub=['sta=~/' char(sta(i)) '/ && orid=~/'];
    for j=1:length(tmpid)-1
        sta_sub=[sta_sub num2str(tmpid(j)) '|'];
    end
    sta_sub=[sta_sub num2str(tmpid(end)) '/'];
    db_tmp=dbsubset(db,sta_sub);
    [orid,orig_time,arr_time]=dbgetv(db_tmp,'orid','origin.time','arrival.time');
    
    if event_from_db
        for j=1:length(orid)
            fprintf(fid_phase,'1 %.2f 1 5\n',arr_time(j)-orig_time(j));
            fprintf(fid_phase,'%d %d 0 1.0 1.0\n',i,orid(j));
        end
    else
        for j=1:length(orid)
            a=find(ev_orid==orid(j));
            fprintf(fid_phase,'1 %.2f 1 5\n',arr_time(j)-(orig_time(j)+ev_time(a)));
            fprintf(fid_phase,'%d %d 0 1.0 1.0\n',i,orid(j));
        end   
    end
end

fclose(fid_phase);

disp('Locating and Inverting')

system(['(/home/menke/raytrace3d/RAYTRACE3D/raytrace3d locate_and_invert ' model_file ' ' model_pars_file ' station.txt ' ...
    locations ' phase.txt inversion_model.model ' num2str(damping) ' ' num2str(rel_weight) ...
    ' > events_updated.txt) >& locate_invert.log']);

dbclose(db);











