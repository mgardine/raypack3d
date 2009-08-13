function [varargout] = ray_plotmodel(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_plotmodel
%
% This function plots the raytrace3d model along the xy, xz, or yz planes.
%
% This function requires the presence of one other function:
%   ray_plotslice.m (for volumetric slice plots)
%       OR
%   ray_plotcontour.m (for contour plots)
%
% It has two uses:
% 
% Usage 1:
% [[X],[Y],[Z],[w]]=ray_plotmodel(r,direction,line,[display],[method],[num_x],[num_y],[num_z])
%
% This usage takes a raytrace3d model and plots along a user-specified
% plane. It can return the gridded data values for later use.
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
% Optional Inputs:
%   display:        The plotting type.  Options are 'slice' (volumetric
%                   slice plot), or 'contour' (contour plot).
%                   (default: 'slice')
%
%   method:         Method for plotting.  Options are 'interp' (interpolate
%                   between nodes), or 'facet' for grid values between nodes.
%                   Note: if method = 'facet', the num_x, num_y, and num_z
%                   values are not used.  Also, 'facet' will not work with
%                   non-parallel grid lines.
%                   (default: 'interp')
%
%   num_x:          The number of x nodes to interpolate in the plot
%                   (default: 50)
%
%   num_y:          The number of y nodes to interpolate in the plot
%                   (default: 50)
%
%   num_z:          The number of z nodes to interpolate in the  plot
%                   (default: 50)
%
%   
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
% ray_plotmodel(X,Y,Z,w,direction,line,[display],[method])
%
% This usage takes the already gridded data output from a previous
% ray_plotmodel run and plots the data along a user-specified plane
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
% Optional Inputs:
%   display:        The plotting type.  Options are 'slice' (volumetric
%                   slice plot), or 'contour' (contour plot).
%                   (default: 'slice')
%
%   method:         Method for plotting.  Options are 'interp' (interpolate
%                   between nodes), or 'facet' for grid values between nodes.
%                   Note: if method = 'facet', the num_x, num_y, and num_z
%                   values are not used.  Also, 'facet' will not work with
%                   non-parallel grid lines.
%                   (default: 'interp')
%
% Author:
% Matt Gardine
% February 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks for the existence of the ray_plotslice function 
if (exist('ray_plotslice') ~= 2)
    error('Error: This function is dependent on ray_plotslice. Please add this function into the path')
end

% Checks for the existence of the ray_plotcontour function 
if (exist('ray_plotcontour') ~= 2)
    error('Error: This function is dependent on ray_plotcontour. Please add this function into the path')
end

if isstruct(varargin{1})
    switch nargin
        case 3
        r = varargin{1};
        num_x=50;
        num_y=50;
        num_z=50;
        direction=varargin{2};
        line=varargin{3};
        display='slice';
        method='interp';
        
        case 4
        if (strcmp(varargin{4},'slice')||strcmp(varargin{4},'contour'))
            r = varargin{1};
            direction=varargin{2};
            line=varargin{3};
            display=varargin{4};
            method='interp';
            num_x=50;
            num_y=50;
            num_z=50;
        elseif (strcmp(varargin{4},'interp')||strcmp(varargin{4},'facet'))
            r = varargin{1};
            direction=varargin{2};
            line=varargin{3};
            method=varargin{4};
            display='slice';
            num_x=50;
            num_y=50;
            num_z=50;
        end
            
        case 5
        r = varargin{1};
        direction=varargin{2};
        line=varargin{3};
        display=varargin{4};
        method=varargin{5};
        num_x=50;
        num_y=50;
        num_z=50;
        
        case 6
        r = varargin{1};
        direction=varargin{2};
        line=varargin{3};
        num_x=varargin{4};
        num_y=varargin{5};
        num_z=varargin{6};
        display='slice';
        method='interp';
        
        case 7
        if (strcmp(varargin{4},'slice')||strcmp(varargin{4},'contour'))
            r = varargin{1};
            direction=varargin{2};
            line=varargin{3};
            display=varargin{4};
            num_x=varargin{5};
            num_y=varargin{6};
            num_z=varargin{7};
            method='interp';
        elseif (strcmp(varargin{4},'interp')||strcmp(varargin{4},'facet'))
            r = varargin{1};
            direction=varargin{2};
            line=varargin{3};
            method=varargin{4};
            num_x=varargin{5};
            num_y=varargin{6};
            num_z=varargin{7};
            display='slice';
        end
        
        case 8
        r = varargin{1};
        direction=varargin{2};
        line=varargin{3};
        display=varargin{4};
        method=varargin{5};
        num_x=varargin{6};
        num_y=varargin{7};
        num_z=varargin{8};
        
        otherwise
        help ray_plotmodel
        return;
    end
    
    if strcmp(method,'facet')
        [X,Y,Z,w]=make_griddata_facet(r);
    else
        [X,Y,Z,w]=make_griddata_interp(r,num_x,num_y,num_z);
    end
    
elseif nargin==6
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    w=varargin{4};
    direction=varargin{5};
    line=varargin{6};
    method = 'interp';
    
elseif nargin==7
    if (strcmp(varargin{7},'slice')||strcmp(varargin{7},'contour'))
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        w=varargin{4};
        direction=varargin{5};
        line=varargin{6};
        display=varargin{7};
        method='interp';
        
    elseif (strcmp(varargin{7},'interp')||strcmp(varargin{7},'facet'))
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        w=varargin{4};
        direction=varargin{5};
        line=varargin{6};
        method=varargin{7};
        display='slice';
    end
    
elseif nargin==8
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    w=varargin{4};
    direction=varargin{5};
    line=varargin{6};
    display=varargin{7};
    method=varargin{8};
    
else
    help ray_plotmodel
    return;
end

if strcmp(display,'slice')
    ray_plotslice(X,Y,Z,w,direction,line,method);

elseif strcmp(display,'contour')
    ray_plotcontour(X,Y,Z,w,direction,line)
end

switch nargout
    case 0
        return
    case 4
        varargout{1}=X;
        varargout{2}=Y;
        varargout{3}=Z;
        varargout{4}=w;
    otherwise
        disp('Error: Invalid number of output arguments');
        return
end



function [X,Y,Z,w]=make_griddata_interp(r,num_x,num_y,num_z)
min_x = min(r.modelxyz(:,1));
max_x = max(r.modelxyz(:,1));
min_y = min(r.modelxyz(:,2));
max_y = max(r.modelxyz(:,2));
min_z = min(r.modelxyz(:,3));
max_z = max(r.modelxyz(:,3));

xi = (min_x:(max_x-min_x)/num_x:max_x);
yi = (min_y:(max_y-min_y)/num_y:max_y);
zi = (min_z:(max_z-min_z)/num_z:max_z);

[X,Y,Z] = meshgrid(xi,yi,zi);
w = griddata3(r.modelxyz(:,1),r.modelxyz(:,2),r.modelxyz(:,3),r.modelxyz(:,4),X,Y,Z,'linear');


function [X,Y,Z,w]=make_griddata_facet(r)
num_x = min(r.modeldims(:,1));
num_y = min(r.modeldims(:,2));
num_z = max(r.modeldims(:,3));

x=zeros(num_x,1);
y=zeros(num_y,1);
z=zeros(num_z,1);

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
end

[X,Y,Z]=meshgrid(x,y,z);
w = griddata3(r.modelxyz(:,1),r.modelxyz(:,2),r.modelxyz(:,3),r.modelxyz(:,4),X,Y,Z,'linear');