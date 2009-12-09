function ray_plotrays(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_plotrays
%
% This function plots the rays to an individual station given by lat, lon
% coordinates, or plots all rays if no lat, lon given
% 
% This function requires the presence of one other function:
%   ray_latlon2xyz.m
%
% Usage: ray_plotrays(r,lat,lon,[color])
%
% Inputs:
%   r:              A matlab raytrace3d structure 
%
%   lat:            The latitude coordinate of the station of interest 
%
%   lon:            The longitude coordinate of the station of interest
%
% Note: Function plots all rays if no lat & lon coordinates given
%
% Optional Input:
%   color:          Which color to plot the rays
%                   (default: blue)
%
% Author:
% Matt Gardine
% November 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
   help ray_plotrays
   return;
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

if isstruct(varargin{1})
    switch nargin
        case 1
            r = varargin{1};
            color = 'b';
            plot_all=1;
            
        case 3
            r = varargin{1};
            lat = varargin{2};
            lon = varargin{3};
            color = 'b';
            plot_all=0;
            
        case 4
            r = varargin{1};
            lat = varargin{2};
            lon = varargin{3};
            color = varargin{4};
            plot_all=0;
            
        otherwise
            help ray_plotrays
            return;
    end

else
    help ray_plotrays
    return;
end

if plot_all
    rays = [r.rayid(:) r.raycoord(:,:)];

else
    [x,y,z]=ray_latlon2xyz(lat,lon,0,ref_lat,ref_lon,projection);

    a = find(r.raycoord(:,1)>x-.005 & r.raycoord(:,1)<x+.005 & r.raycoord(:,2)>y-.005 & r.raycoord(:,2)<y+.005);

    k=1;
    for i=1:length(a)
        for j=a(i):-1:1
            if r.rayid(j)==r.rayid(a(i))
                b(k)=j;
                k=k+1;
            else
                break
            end
        end
    end

    rays = [r.rayid(b) r.raycoord(b,:)];
end

rayid = rays(1,1);
temprays = rays(1,:);
count=0;

%figure
hold on

for i=2:length(rays)
    if rays(i,1)==rayid && rays(i,4)+10>rays(i-1,4)
        temprays=[temprays; rays(i,:)];
    else
        plot3(temprays(end,2),temprays(end,3),temprays(end,4),'ko','MarkerFaceColor','k','MarkerSize',6)
        plot3(temprays(2:end,2),temprays(2:end,3),temprays(2:end,4),[color '-'],'LineWidth',1)
        plot3(temprays(1,2),temprays(1,3),temprays(1,4),'kv','MarkerFaceColor','k','MarkerSize',8)
        temprays=rays(i,:);
        rayid=rays(i,1);
        count = count+1;
    end
end
disp([num2str(count) ' rays plotted'])