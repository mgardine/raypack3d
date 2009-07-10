function r = ray_loaddata(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_loaddata
%
% This function loads data from raytrace3d traveltime file, ray file, model
% file, and model hits file into a matlab raytrace3d structure
% 
% Usage: r = ray_loaddata(model_file,[tt_file],[location_file],[ray_file],[model_hits_file])
%
% Required Input:
%   model_file:     The 3D model file used in any raytrace3d routine
%
% Optional Inputs:
%   tt_file:        The travel time file created by running 
%                   ray_shoot_predicted
%
%   location_file:  The hypocentral origin file created by running
%                   ray_locate
%
%   ray_file:       The ray file created by running ray_shoot_predicted
%
%   model_hits_file: The 3D model file created by the ray_invert_db or
%                    ray_invert_ttfile, corresponding to the number of rays
%                    passing through each grid node
%
% Output:
%   r:              A matlab raytrace3d structure with the following fields:
%
%                   hypocenter: [X, Y, Z, origin time] of earthquakes
%                   station: [X, Y, Z, arrival time] of stations
%                   traveltime: Arrival time - origin time
%                   residual: Observed arrival time - predicted arrival time
%                   rayid: Ray identification number
%                   raycoord: [X, Y, Z] coordinates of the ray
%                   tetra: Tetrahedron number that ray travels through
%                   modeldims: [Num_X Num_Y Num_Z] number of nodes in model
%                   modelijk: [node i, node j, node k, velocity, interface, hits]
%                   modelxyz: [node x, node y, node z, veloctiy, interface, hits]
%
% Author:
% Matt Gardine
% February 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializes optional fields to blank values
x1=[];
y1=[];
z1=[];
Ta=[];
TT=[];
dT=[];
idr=[];
Xr=[];
Yr=[];
Zr=[];
tetrar=[];

switch nargin
% Checks if only the required model file is given
case 1
    [numi,numj,numk]=textread(varargin{1},'%d %d %d',1);
    model_file=varargin{1};
    [i j k x y z v int]=textread(model_file,'%d %d %d %f %f %f %f %d','headerlines',1);
    x = round(x./0.01).*0.01;
    y = round(y./0.01).*0.01;
    z = round(z./0.01).*0.01;
    v = round(v./0.01).*0.01;
    for temp=1:length(i)
        hits(temp)=0;
    end
    hits=hits';
    
% Checks if the required model file plus one other file is given
case 2
    [numi,numj,numk]=textread(varargin{1},'%d %d %d',1);
    model_file=varargin{1};
    [i j k x y z v int]=textread(model_file,'%d %d %d %f %f %f %f %d','headerlines',1);
    x = round(x./0.01).*0.01;
    y = round(y./0.01).*0.01;
    z = round(z./0.01).*0.01;
    v = round(v./0.01).*0.01;
    
    for temp=1:length(i)
        hits(temp)=0;
    end
    hits=hits';
    
    [s1,s2] = textread(varargin{2},'%s %s',1);
    
    % Checks if the second argument is a travel-time file
    if strcmp(s2,'x0')
        tt_file=varargin{2};
        [id x0 y0 z0 x1 y1 z1 T0 Ta TT T dT ratio obs_order pre_order amplitude ntts]=textread(tt_file,'%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d','headerlines',1);
     
    % Checks if the second argument is a location file
    elseif strcmp(s2,'xs')
        location_file=varargin{2};
        [orid xs ys zs t0 sigx sigy sigz sigt ndat iter rms]=textread(location_file,'%d %f %f %f %f %f %f %f %f %d %d %f','headerlines',1);
        
    % Checks if the second argument is a ray file
    elseif strcmp(s2,'x')
        ray_file=varargin{2};
        [idr Xr Yr Zr tetrar ]=textread(ray_file,'%d %f %f %f %d','headerlines',1);
        
    % Otherwise, assumes that the last file is a model.hits file
    else
        model_hits_file=varargin{2};
        [junk1 junk2 junk3 junk4 junk5 junk6 hits junk7]=textread(model_hits_file,'%d %d %d %f %f %f %d %d','headerlines',1);
    end

% Checks if the require model file, plus two other files are given    
case 3
    [numi,numj,numk]=textread(varargin{1},'%d %d %d',1);
    model_file=varargin{1};
    [i j k x y z v int]=textread(model_file,'%d %d %d %f %f %f %f %d','headerlines',1);
    x = round(x./0.01).*0.01;
    y = round(y./0.01).*0.01;
    z = round(z./0.01).*0.01;
    v = round(v./0.01).*0.01;
    
    for temp=1:length(i)
        hits(temp)=0;
    end
    hits=hits';
    
    [s1,s2] = textread(varargin{2},'%s %s',1);
    [s3,s4] = textread(varargin{3},'%s %s',1);
    
    % Checks what file type the second input is
    % Checks if second input is a travel time file
    if strcmp(s2,'x0')
        tt_file=varargin{2};
        [id x0 y0 z0 x1 y1 z1 T0 Ta TT T dT ratio obs_order pre_order amplitude ntts]=textread(tt_file,'%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d','headerlines',1);
    
    % Checks if the second input is a location file
    elseif strcmp(s2,'xs')
        location_file=varargin{2};
        [orid xs ys zs t0 sigx sigy sigz sigt ndat iter rms]=textread(location_file,'%d %f %f %f %f %f %f %f %f %d %d %f','headerlines',1);
           
    % Checks if the second input is a ray file
    elseif strcmp(s2,'x')
        ray_file=varargin{2};
        [idr Xr Yr Zr tetrar ]=textread(ray_file,'%d %f %f %f %d','headerlines',1);
        
    % Otherwise, assumes that the second input is a model.hits file
    else
        model_hits_file=varargin{2};
        [junk1 junk2 junk3 junk4 junk5 junk6 hits junk7]=textread(model_hits_file,'%d %d %d %f %f %f %d %d','headerlines',1);
    end
    
    
    % Checks what file type the third input is
    % Checks if third input is a travel time file
    if strcmp(s4,'x0')
        tt_file=varargin{3};
        [id x0 y0 z0 x1 y1 z1 T0 Ta TT T dT ratio obs_order pre_order amplitude ntts]=textread(tt_file,'%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d','headerlines',1);
    
    % Checks if the third input is a location file
    elseif strcmp(s4,'xs')
        location_file=varargin{3};
        [orid xs ys zs t0 sigx sigy sigz sigt ndat iter rms]=textread(location_file,'%d %f %f %f %f %f %f %f %f %d %d %f','headerlines',1);
              
    % Checks if the third input is a ray file
    elseif strcmp(s4,'x')
        ray_file=varargin{3};
        [idr Xr Yr Zr tetrar ]=textread(ray_file,'%d %f %f %f %d','headerlines',1);
    
    % Otherwise, assumes that the third input is a model.hits file
    else
        model_hits_file=varargin{3};
        [junk1 junk2 junk3 junk4 junk5 junk6 hits junk7]=textread(model_hits_file,'%d %d %d %f %f %f %d %d','headerlines',1);
    end
 
    
% Checks if four files are given
case 4
    [numi,numj,numk]=textread(varargin{1},'%d %d %d',1);
    model_file=varargin{1};
    [i j k x y z v int]=textread(model_file,'%d %d %d %f %f %f %f %d','headerlines',1);
    x = round(x./0.01).*0.01;
    y = round(y./0.01).*0.01;
    z = round(z./0.01).*0.01;
    v = round(v./0.01).*0.01;
    
    [s1,s2] = textread(varargin{2},'%s %s',1);
    [s3,s4] = textread(varargin{3},'%s %s',1);
    [s5,s6] = textread(varargin{4},'%s %s',1);
    
    % Checks what file type the second input is
    % Checks if second input is a travel time file
    if strcmp(s2,'x0')
        tt_file=varargin{2};
        [id x0 y0 z0 x1 y1 z1 T0 Ta TT T dT ratio obs_order pre_order amplitude ntts]=textread(tt_file,'%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d','headerlines',1);
    
    % Checks if the second input is a location file
    elseif strcmp(s2,'xs')
        location_file=varargin{2};
        [orid xs ys zs t0 sigx sigy sigz sigt ndat iter rms]=textread(location_file,'%d %f %f %f %f %f %f %f %f %d %d %f','headerlines',1);
            
    % Checks if the second input is a ray file
    elseif strcmp(s2,'x')
        ray_file=varargin{2};
        [idr Xr Yr Zr tetrar ]=textread(ray_file,'%d %f %f %f %d','headerlines',1);
        
    % Otherwise, assumes that the second input is a model.hits file
    else
        model_hits_file=varargin{2};
        [junk1 junk2 junk3 junk4 junk5 junk6 hits junk7]=textread(model_hits_file,'%d %d %d %f %f %f %d %d','headerlines',1);
    end
    
    
    % Checks what file type the third input is
    % Checks if third input is a travel time file
    if strcmp(s4,'x0')
        tt_file=varargin{3};
        [id x0 y0 z0 x1 y1 z1 T0 Ta TT T dT ratio obs_order pre_order amplitude ntts]=textread(tt_file,'%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d','headerlines',1);
    
    % Checks if the third input is a location file
    elseif strcmp(s4,'xs')
        location_file=varargin{3};
        [orid xs ys zs t0 sigx sigy sigz sigt ndat iter rms]=textread(location_file,'%d %f %f %f %f %f %f %f %f %d %d %f','headerlines',1);
           
    % Checks if the third input is a ray file
    elseif strcmp(s4,'x')
        ray_file=varargin{3};
        [idr Xr Yr Zr tetrar ]=textread(ray_file,'%d %f %f %f %d','headerlines',1);
    
    % Otherwise, assumes that the third input is a model.hits file
    else
        model_hits_file=varargin{3};
        [junk1 junk2 junk3 junk4 junk5 junk6 hits junk7]=textread(model_hits_file,'%d %d %d %f %f %f %d %d','headerlines',1);
    end   
    
    % Checks what file type the fourth input is
    % Checks if fourth input is a travel time file
    if strcmp(s6,'x0')
        tt_file=varargin{4};
        [id x0 y0 z0 x1 y1 z1 T0 Ta TT T dT ratio obs_order pre_order amplitude ntts]=textread(tt_file,'%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d','headerlines',1);
        
    % Checks if the third input is a location file
    elseif strcmp(s6,'xs')
        location_file=varargin{4};
        [orid xs ys zs t0 sigx sigy sigz sigt ndat iter rms]=textread(location_file,'%d %f %f %f %f %f %f %f %f %d %d %f','headerlines',1);
     
    % Checks if the fourth input is a ray file
    elseif strcmp(s6,'x')
        ray_file=varargin{4};
        [idr Xr Yr Zr tetrar ]=textread(ray_file,'%d %f %f %f %d','headerlines',1);
    
    % Otherwise, assumes that the fourth input is a model.hits file
    else
        model_hits_file=varargin{4};
        [junk1 junk2 junk3 junk4 junk5 junk6 hits junk7]=textread(model_hits_file,'%d %d %d %f %f %f %d %d','headerlines',1);
    end 
    
% Checks if all five files are given
case 5
    [numi,numj,numk]=textread(varargin{1},'%d %d %d',1);
    model_file=varargin{1};
    [i j k x y z v int]=textread(model_file,'%d %d %d %f %f %f %f %d','headerlines',1);
    x = round(x./0.01).*0.01;
    y = round(y./0.01).*0.01;
    z = round(z./0.01).*0.01;
    v = round(v./0.01).*0.01;
    
    [s1,s2] = textread(varargin{2},'%s %s',1);
    [s3,s4] = textread(varargin{3},'%s %s',1);
    [s5,s6] = textread(varargin{4},'%s %s',1);
    [s7,s8] = textread(varargin{5},'%s %s',1);
    
    % Checks what file type the second input is
    % Checks if second input is a travel time file
    if strcmp(s2,'x0')
        tt_file=varargin{2};
        [id x0 y0 z0 x1 y1 z1 T0 Ta TT T dT ratio obs_order pre_order amplitude ntts]=textread(tt_file,'%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d','headerlines',1);
    
    % Checks if the second input is a location file
    elseif strcmp(s2,'xs')
        location_file=varargin{2};
        [orid xs ys zs t0 sigx sigy sigz sigt ndat iter rms]=textread(location_file,'%d %f %f %f %f %f %f %f %f %d %d %f','headerlines',1);
            
    % Checks if the second input is a ray file
    elseif strcmp(s2,'x')
        ray_file=varargin{2};
        [idr Xr Yr Zr tetrar ]=textread(ray_file,'%d %f %f %f %d','headerlines',1);
        
    % Otherwise, assumes that the second input is a model.hits file
    else
        model_hits_file=varargin{2};
        [junk1 junk2 junk3 junk4 junk5 junk6 hits junk7]=textread(model_hits_file,'%d %d %d %f %f %f %d %d','headerlines',1);
    end
    
    
    % Checks what file type the third input is
    % Checks if third input is a travel time file
    if strcmp(s4,'x0')
        tt_file=varargin{3};
        [id x0 y0 z0 x1 y1 z1 T0 Ta TT T dT ratio obs_order pre_order amplitude ntts]=textread(tt_file,'%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d','headerlines',1);
    
    % Checks if the third input is a location file
    elseif strcmp(s4,'xs')
        location_file=varargin{3};
        [orid xs ys zs t0 sigx sigy sigz sigt ndat iter rms]=textread(location_file,'%d %f %f %f %f %f %f %f %f %d %d %f','headerlines',1);
           
    % Checks if the third input is a ray file
    elseif strcmp(s4,'x')
        ray_file=varargin{3};
        [idr Xr Yr Zr tetrar ]=textread(ray_file,'%d %f %f %f %d','headerlines',1);
    
    % Otherwise, assumes that the third input is a model.hits file
    else
        model_hits_file=varargin{3};
        [junk1 junk2 junk3 junk4 junk5 junk6 hits junk7]=textread(model_hits_file,'%d %d %d %f %f %f %d %d','headerlines',1);
    end   
    
    % Checks what file type the fourth input is
    % Checks if fourth input is a travel time file
    if strcmp(s6,'x0')
        tt_file=varargin{4};
        [id x0 y0 z0 x1 y1 z1 T0 Ta TT T dT ratio obs_order pre_order amplitude ntts]=textread(tt_file,'%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d','headerlines',1);
        
    % Checks if the fourth input is a location file
    elseif strcmp(s6,'xs')
        location_file=varargin{4};
        [orid xs ys zs t0 sigx sigy sigz sigt ndat iter rms]=textread(location_file,'%d %f %f %f %f %f %f %f %f %d %d %f','headerlines',1);
     
    % Checks if the fourth input is a ray file
    elseif strcmp(s6,'x')
        ray_file=varargin{4};
        [idr Xr Yr Zr tetrar ]=textread(ray_file,'%d %f %f %f %d','headerlines',1);
    
    % Otherwise, assumes that the fourth input is a model.hits file
    else
        model_hits_file=varargin{4};
        [junk1 junk2 junk3 junk4 junk5 junk6 hits junk7]=textread(model_hits_file,'%d %d %d %f %f %f %d %d','headerlines',1);
    end     
    
    % Checks what file type the fifth input is
    % Checks if fifth input is a travel time file
    if strcmp(s8,'x0')
        tt_file=varargin{5};
        [id x0 y0 z0 x1 y1 z1 T0 Ta TT T dT ratio obs_order pre_order amplitude ntts]=textread(tt_file,'%d %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %d','headerlines',1);
        
    % Checks if the fifth input is a location file
    elseif strcmp(s8,'xs')
        location_file=varargin{5};
        [orid xs ys zs t0 sigx sigy sigz sigt ndat iter rms]=textread(location_file,'%d %f %f %f %f %f %f %f %f %d %d %f','headerlines',1);
     
    % Checks if the fifth input is a ray file
    elseif strcmp(s8,'x')
        ray_file=varargin{5};
        [idr Xr Yr Zr tetrar ]=textread(ray_file,'%d %f %f %f %d','headerlines',1);
    
    % Otherwise, assumes that the fifth input is a model.hits file
    else
        model_hits_file=varargin{5};
        [junk1 junk2 junk3 junk4 junk5 junk6 hits junk7]=textread(model_hits_file,'%d %d %d %f %f %f %d %d','headerlines',1);
    end 

end

if (exist('xs','var')==1 && exist('x0','var')==1)
    % If both location_file and tt_file exist, use the origins from the tt_file
    display('Both location_file and tt_file exist, using origins from tt_file')
elseif exist('xs','var')==1
    x0=xs;
    y0=ys;
    z0=zs;
    T0=t0;
else
    x0=[];
    y0=[];
    z0=[];
    T0=[];
end


% Creates the raytrace structure
r = struct('hypocenter',[],'station',[],'traveltime',[],'residual',[],'rayid',[],'raycoord',[],'tetra',[],'modeldims',[],'modelijk',[],'modelxyz',[]);

r.hypocenter = [x0 y0 z0 T0];
r.station = [x1 y1 z1 Ta];
r.traveltime = TT;
r.residual = dT;
r.rayid = idr;
r.raycoord = [Xr Yr Zr];
r.tetra = tetrar;
r.modeldims = [numi numj numk];
r.modelijk = [i j k v int hits];
r.modelxyz = [x y z v int hits];

% Sorts the model rows into a common format
[tmp index]=sortrows([r.modelijk(:,1) r.modelijk(:,2) r.modelijk(:,3)],[3 2 1]);
r.modelijk(:,:)=r.modelijk(index,:);
r.modelxyz(:,:)=r.modelxyz(index,:);


