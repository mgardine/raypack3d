function [varargout] = ray_plothits(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_plothits
%
% This function plots the number of ray hits in each grid node in a 
% raytrace3d model along the xy, xz, or yz planes.  To maintain a
% meaningful color scale, the maximum number of ray hits per grid node is
% set to be 100.
%
% This function requires the presence of one other function:
%   ray_plotslice.m
%
% It has two uses:
% 
% Usage 1:
% [[X],[Y],[Z],[w]]=ray_plothits(r,direction,line)
%
% This usage takes a raytrace3d structure and plots along the number of ray
% hits in each grid node along a user-specified plane.  It can return the
% gridded data values for later use.
%
% Required Inputs:
%   r:              A matlab raytrace3d structure 
%
%   direction:      Determines which plane to plot. Valid options are 
%                   'xy', 'xz', or 'yz'
%
%   line:           The value along the third dimension to create the
%                   volume slice through   
%
%
% Optional Outputs:
%   X:              X grid nodes output from meshgrid
%
%   Y:              Y grid nodes, output from meshgrid
%
%   Z:              Z grid nodes, output from meshgrid
%
%   w:              Data at each grid node, output from griddata3
%
%
% Usage 2: 
% ray_plothits(X,Y,Z,w,direction,line)
%
% This usage takes the already gridded data output from a previous
% ray_plothits run and plots the data along a user-specified plane
%
% Inputs:
%   X:              X grid nodes output from meshgrid
%
%   Y:              Y grid nodes, output from meshgrid
%
%   Z:              Z grid nodes, output from meshgrid
%
%   w:              Data at each grid node, output from griddata3
%
%   direction:      Determines which plane to plot. Valid options are 
%                   'xy', 'xz', or 'yz'
%
%   line:           The value along the third dimension to create the
%                   volume slice through
%
%
% Author:
% Matt Gardine
% May 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks for the existence of the ray_plotslice function 
if (exist('ray_plotslice') ~= 2)
    error('Error: This function is dependent on ray_plotslice. Please add this function into the path')
end

if isstruct(varargin{1})
    if nargin==3
        r = varargin{1};
        direction=varargin{2};
        line=varargin{3};
    else
        disp('Error: Invalid input parameters');
        return
    end
    
    [X,Y,Z,w]=make_griddata_facet(r);
    
elseif nargin==6
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    w=varargin{4};
    direction=varargin{5};
    line=varargin{6};
    
else
    disp('Error: Input must be a raytrace3d structure or gridded data');
    return
end

ray_plotslice(X,Y,Z,w,direction,line,'facet');

if nargout==4
    varargout{1}=X;
    varargout{2}=Y;
    varargout{3}=Z;
    varargout{4}=w;
end


function [X,Y,Z,w]=make_griddata_facet(r)
num_x = min(r.modeldims(:,1));
num_y = min(r.modeldims(:,2));
num_z = max(r.modeldims(:,3));

j=1;
k=1;
m=1;
n=1;

for i=1:length(r.modelxyz(:,6))    
    if i==j && i<(num_x*num_z)
        x(k)=r.modelxyz(i,1);
        j=(num_z*k)+1;
        k=k+1;
    end
    if i==m
        y(n)=r.modelxyz(i,2);
        m=(num_x*num_z*n)+1;
        n=n+1;
    end
    if (i-1)<num_z
        z(i)=r.modelxyz(i,3);
    end
    if r.modelxyz(i,6)>100
        r.modelxyz(i,6)=100;
    end
    
end

[X,Y,Z]=meshgrid(x,y,z);
w = griddata3(r.modelxyz(:,1),r.modelxyz(:,2),r.modelxyz(:,3),r.modelxyz(:,6),X,Y,Z,'linear');

