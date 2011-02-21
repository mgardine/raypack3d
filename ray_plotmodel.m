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
% [[X],[Y],[Z],[w],[handle]]=ray_plotmodel(r,direction,line,[display],[interval],[method],[num_x],[num_y],[num_z],[datatype])
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
%   interval:       Sets a user-defined contour interval. This input is
%                   only used when 'display' is set to 'contour'. Note
%                   that, if used, this input must immediately follow the
%                   declaration of 'contour'.
%                   (default: 0.2 km/s)
%
%   method:         Method for plotting.  Options are 'interp' (interpolate
%                   between nodes), or 'facet' for grid values between nodes.
%                   Note: if method = 'facet', the num_x, num_y, and num_z
%                   values are not used.  Also, 'facet' will not work with
%                   non-parallel grid lines.
%                   (default: 'interp')
%
%   num_x:          The number of x nodes to interpolate in the plot. Note
%                   that, if used, this input must immediately follow the
%                   declaration of 'interp' and be followed by num_y and num_z.
%                   (default: 50)
%
%   num_y:          The number of y nodes to interpolate in the plot. Note
%                   that this input must immediately follow the value given
%                   for num_x, and be followed by num_z. 
%                   (default: 50)
%
%   num_z:          The number of z nodes to interpolate in the plot. Note
%                   that this input must immediately follow the value given
%                   by num_y.
%                   (default: 50)
%
%   datatype:       Datatype from raypack3d structure to plot. Options are
%                   'velocity' to plot velocities (from r.modelxyz(:,4), 
%                   values filled from ray_invert_db), or 'PS' to plot
%                   either Q or P:S ratio (from r.modelxyz(:,7), values
%                   filled from ray_invert_q or ray_invert_ps).
%                   (default: 'velocity')
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
%   handle:         The graphics handle for the given plot, useful for 
%                   low-level editing of the figure.
%
% Usage 2: 
% [handle]=ray_plotmodel(X,Y,Z,w,direction,line,[display],[interval],[method])
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
%   interval:       Sets a user-defined contour interval. This input is
%                   only used when 'display' is set to 'contour'.  Note
%                   that, if used, this input must immediately follow the
%                   declaration of 'contour'.
%                   (default: 0.2 km/s)  
%
%   method:         Method for plotting.  Options are 'interp' (interpolate
%                   between nodes), or 'facet' for grid values between nodes.
%                   Note: if method = 'facet', the num_x, num_y, and num_z
%                   values are not used.  Also, 'facet' will not work with
%                   non-parallel grid lines.
%                   (default: 'interp')
%
%
% Optional Output:
%   handle:         The graphics handle for the given plot, useful for 
%                   low-level editing of the figure .
%
% Author:
% Matt Gardine
% February 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
   help ray_plotmodel
   return;
end

% Checks for the existence of the ray_plotslice function 
if (exist('ray_plotslice') ~= 2)
    error('Error: This function is dependent on ray_plotslice. Please add this function into the path')
end

% Checks for the existence of the ray_plotcontour function 
if (exist('ray_plotcontour') ~= 2)
    error('Error: This function is dependent on ray_plotcontour. Please add this function into the path')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if a raypack3d structure is the first argument
if isstruct(varargin{1})
    switch nargin
        % Checks for the minimum input arguments
        case 3
        r = varargin{1};
        num_x=50;
        num_y=50;
        num_z=50;
        direction=varargin{2};
        line=varargin{3};
        display='slice';
        method='interp';
        datatype='velocity';
        
        
        % Checks if one optional arguments is given
        case 4
        % Checks if the first optional argument is for DISPLAY 
        if (strcmp(varargin{4},'slice')||strcmp(varargin{4},'contour'))
            r = varargin{1};
            direction=varargin{2};
            line=varargin{3};
            display=varargin{4};
            if strcmp(display,'contour')
                interval=0.2;
            end
            method='interp';
            datatype='velocity';
            num_x=50;
            num_y=50;
            num_z=50;
        % Checks if the first optional argument is for METHOD
        elseif (strcmp(varargin{4},'interp')||strcmp(varargin{4},'facet'))
            r = varargin{1};
            direction=varargin{2};
            line=varargin{3};
            method=varargin{4};
            display='slice';
            datatype='velocity';
            num_x=50;
            num_y=50;
            num_z=50;
        % Checks if the first optional argument is for DATATYPE
        elseif (strcmp(varargin{4},'velocity')||strcmp(varargin{4},'ps'))
            r = varargin{1};
            direction=varargin{2};
            line=varargin{3};
            datatype=varargin{4};
            display='slice';
            method='interp';
            num_x=50;
            num_y=50;
            num_z=50;
        % Otherwise, show help file and quit routine
        else
            help ray_plotmodel
            return;
        end
        
        
        % Checks if two optional arguments are given
        case 5
        r = varargin{1};
        direction=varargin{2};
        line=varargin{3};
        
        % Checks if the first optional argument is for DISPLAY
        if strcmp(varargin{4},'slice')||strcmp(varargin{4},'contour')
            display=varargin{4};
            % Checks if the second optional argument is for METHOD
            if strcmp(varargin{5},'interp')||strcmp(varargin{5},'facet')
                method=varargin{5};
                interval=0.2;
                datatype='velocity';
            % Checks if the second optional argument is for DATATYPE
            elseif strcmp(varargin{5},'velocity')||strcmp(varargin{5},'ps')
                datatype=varargin{5};
                method='interp';
                interval=0.2;
            % Checks if the second optional argument is for INTERVAL
            elseif strcmp(display,'contour')&&isnumeric(varargin{5})
                interval=varargin{5};
                method='interp';
                datatype='velocity';
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
            
        % Checks if the first optional argument is for DATATYPE    
        elseif strcmp(varargin{4},'velocity')||strcmp(varargin{4},'ps')
            datatype=varargin{4};
            % Checks if the second optional argument is for METHOD
            if strcmp(varargin{5},'interp')||strcmp(varargin{5},'facet')
                method=varargin{5};
                interval=0.2;
                display='slice';
            % Checks if the second optional argument is for DISPLAY
            elseif strcmp(varargin{5},'slice')||strcmp(varargin{5},'contour')
                display=varargin{5};
                interval=0.2;
                method='interp';
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
            
        % Checks if the first optional argument is for METHOD
        elseif strcmp(varargin{4},'interp')||strcmp(varargin{4},'facet')
            method=varargin{4};
            % Checks if second optional argument is for DISPLAY
            if strcmp(varargin{5},'slice')||strcmp(varargin{5},'contour')
                display=varargin{5};
                interval=0.2;
                datatype='velocity';
            % Checks if second optional argument is for DATATYPE
            elseif strcmp(varargin{5},'velocity')||strcmp(varargin{5},'ps')
                datatype=varargin{5};
                display='slice';
                interval=0.2;
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
        end
            
        num_x=50;
        num_y=50;
        num_z=50;
        
        
        % Checks if three optional arguments are given
        case 6
        r = varargin{1};
        direction=varargin{2};
        line=varargin{3};
        
        % Checks if first optional argument is for DISPLAY
        if strcmp(varargin{4},'slice')||strcmp(varargin{4},'contour')
            display=varargin{4};
            % Checks if second optional argument is for METHOD
            if strcmp(varargin{5},'interp')||strcmp(varargin{5},'facet')
                method=varargin{5};
                % Checks if third optional argument is for DATATYPE
                if strcmp(varargin{6},'velocity')||strcmp(varargin{6},'ps')
                    datatype=varargin{6};
                    interval=0.2;
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
                
            % Checks if second optional argument is for DATATYPE
            elseif strcmp(varargin{5},'velocity')||strcmp(varargin{5},'ps')
                datatype=varargin{5};
                % Checks if third optional argument is for METHOD
                if strcmp(varargin{6},'interp')||strcmp(varargin{6},'facet')
                    method=varargin{6};
                    interval=0.2;
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
                
            % Checks if second optional argument is for INTERVAL
            elseif strcmp(display,'contour')&&isnumeric(varargin{5});
                interval=varargin{5};
                % Checks if third optional argument is for METHOD
                if strcmp(varargin{6},'interp')||strcmp(varargin{6},'facet')
                    method=varargin{6};
                    datatype='velocity';
                % Checks if third optional argument is for DATATYPE
                elseif strcmp(varargin{6},'velocity')||strcmp(varargin{6},'ps')
                    datatype=varargin{6};
                    method='interp';
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
                
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
            num_x=50;
            num_y=50;
            num_z=50;
           
        % Checks if first optional argument is for METHOD
        elseif strcmp(varargin{4},'interp')||strcmp(varargin{4},'facet')
            method=varargin{4};
            % Checks if second optional argument is for DISPLAY
            if strcmp(varargin{5},'slice')||strcmp(varargin{5},'contour')
                display=varargin{5};
                % Checks if third optional argument is for DATATYPE
                if strcmp(varargin{6},'velocity')||strcmp(varargin{6},'ps')
                    datatype=varargin{6};
                    interval=0.2;
                % Checks if third optional argument is for INTERVAL
                elseif strcmp(display,'contour')&&isnumeric(varargin{6})
                    interval=varargin{6};
                    datatype='velocity';
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Checks if second optional argument is for DATATYPE
            elseif strcmp(varargin{5},'velocity')||strcmp(varargin{5},'ps')
                datatype=varargin{5};
                % Checks if third optional argument is for DISPLAY
                if strcmp(varargin{6},'slice')||strcmp(varargin{6},'contour')
                    display=varargin{6};
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
                interval=0.2;
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
            num_x=50;
            num_y=50;
            num_z=50;
        
        % Checks if first optional argument is for DATATYPE
        elseif strcmp(varargin{4},'velocity')||strcmp(varargin{4},'ps')
            datatype=varargin{4};
            % Checks if second optional argument is for DISPLAY
            if strcmp(varargin{5},'slice')||strcmp(varargin{5},'contour')
                display=varargin{5};
                % Checks if third optional argument is for METHOD
                if strcmp(varargin{6},'interp')||strcmp(varargin{6},'facet')
                    method=varargin{6};
                    interval=0.2;
                % Checks if the third optional argument is for INTERVAL
                elseif strcmp(display,'contour')&&isnumeric(varargin{6})
                    interval=varargin{6};
                    method='interp';
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Checks if second optional argument is for METHOD
            elseif strcmp(varargin{5},'interp')||strcmp(varargin{5},'facet')
                method=varargin{5};
                % Checks if third optional argument is for DISPLAY
                if strcmp(varargin{6},'slice')||strcmp(varargin{6},'contour')
                    display=varargin{6};
                    interval=0.2;
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
            num_x=50;
            num_y=50;
            num_z=50;
                
        % Otherwise, show help file and quit routine
        else
            help ray_plotmodel
            return;
        end
        

        % Checks if four optional arguments are given
        case 7
            
        % Checks if first optional argument is for DISPLAY
        if (strcmp(varargin{4},'slice')||strcmp(varargin{4},'contour'))
            r = varargin{1};
            direction=varargin{2};
            line=varargin{3};
            display=varargin{4};
            % Checks if second optional argument is for INTERVAL
            if strcmp(display,'contour')&&isnumeric(varargin{5})
                interval=varargin{5};
                % Checks if third optional argument is for METHOD
                if strcmp(varargin{6},'interp')||strcmp(varargin{6},'facet')
                    method=varargin{6};
                    % Checks if fourth optional argument is for DATATYPE
                    if strcmp(varargin{7},'velocity')||strcmp(varargin{7},'ps')
                        datatype=varargin{7};
                    % Otherwise, show help file and quit routine
                    else
                        help ray_plotmodel
                        return;
                    end
                % Checks if third optional argument is for DATATYPE
                elseif strcmp(varargin{6},'velocity')||strcmp(varargin{6},'ps')
                    datatype=varargin{6};
                    % Checks if fourth optional argument is for METHOD
                    if strcmp(varargin{7},'interp')||strcmp(varargin{7},'facet')
                        method=varargin{7};
                    % Otherwise, show help file and quit routine  
                    else
                        help ray_plotmodel
                        return;
                    end
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Otherwise, show help file and quit routine    
            else
                help ray_plotmodel
                return;
            end
            num_x=50;
            num_y=50;
            num_z=50;
            
        % Checks if first optional argument is for METHOD        
        elseif strcmp(varargin{4},'interp')||strcmp(varargin{4},'facet')
            method=varargin{4};
            % Checks if second, third, and fourth optional arguments are
            % NUM_X, NUM_Y, and NUM_Z
            if strcmp(method,'interp')&&isnumeric(varargin{5})&&isnumeric(varargin{6})&&isnumeric(varargin{7})
                num_x=varargin{5};
                num_y=varargin{6};
                num_z=varargin{7};
                display='slice';
                datatype='velocity';
            % Checks if second optional argument is DISPLAY
            elseif strcmp(varargin{5},'contour')
                display=varargin{5};
                % Checks if third optional argument is INTERVAL
                if isnumeric(varargin{6})
                    interval=varargin{6};
                    datatype=varargin{7};
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Checks if second optional argument is DATATYPE
            elseif strcmp(varargin{5},'velocity')||strcmp(varargin{5},'ps')
                datatype=varargin{5};
                % Checks if third optional argument is DISPLAY
                if strcmp(varargin{6},'contour')
                    display=varargin{6};
                    % Checks if fourth optional argument is INTERVAL
                    if isnumeric(varargin{7})
                        interval=varargin{7};
                    % Otherwise, show help file and quit routine
                    else
                        help ray_plotmodel
                        return;
                    end
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
            num_x=50;
            num_y=50;
            num_z=50;
            
        % Checks if first optional argument is DATATYPE
        elseif strcmp(varargin{4},'velocity')||strcmp(varargin{4},'ps')
            datatype=varargin{4};
            % Checks if second optional argument is METHOD
            if strcmp(varargin{5},'interp')||strcmp(varargin{5},'facet')
                method=varargin{5};
                % Checks if third and fourth optional arguments are DISPLAY
                % and INTERVAL
                if strcmp(varargin{6},'contour')&&isnumeric(varargin{7})
                    display=varargin{6};
                    interval=varargin{7};
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Checks if second, third, and fourth optional arguments are
            % DISPLAY, INTERVAL, and METHOD
            elseif strcmp(varargin{5},'contour')&&isnumeric(varargin{6})&&(strcmp(varargin{7},'interp')||strcmp(varargin{7},'facet'))
                display=varargin{5};
                interval=varargin{6};
                method=varargin{7};
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
            num_x=50;
            num_y=50;
            num_z=50;
            
        % Otherwise, show help file and quit routine
        else
            help ray_plotmodel
            return;
        end
 
             
        
        % Checks if five optional arguments are given
        case 8
        r = varargin{1};
        direction=varargin{2};
        line=varargin{3};
        
        % Checks if the first optional argument is METHOD
        if strcmp(varargin{4},'interp')
            % Checks that the second, third, and fourth optional arguments
            % are NUM_X,NUM_Y,and NUM_Z
            if isnumeric(varargin{5})&&isnumeric(varargin{6})&&isnumeric(varargin{7})
                num_x=varargin{5};
                num_y=varargin{6};
                num_z=varargin{7};
                % Checks if fifth optional argument is DISPLAY
                if strcmp(varargin{8},'slice')||strcmp(varargin{8},'contour')
                    display=varargin{8};
                    interval=0.2;
                    datatype='velocity';
                % Checks if fifth optional argument is DATATYPE
                elseif strcmp(varargin{8},'velocity')||strcmp(varargin{8},'ps')
                    datatype=varargin{8};
                    display='slice';
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
            
        % Checks if the first optional argument is DISPLAY   
        elseif strcmp(varargin{4},'slice')||strcmp(varargin{4},'contour')
            display=varargin{4};
            % Checks that the second, third, fourth, and fifth optional
            % arguments are METHOD, NUM_X, NUM_Y, and NUM_Z
            if strcmp(varargin{5},'interp')&&isnumeric(varargin{6})&&isnumeric(varargin{7})&&isnumeric(varargin{8})
                method=varargin{5};
                num_x=varargin{6};
                num_y=varargin{7};
                num_z=varargin{8};
                datatype='velocity';
                interval=0.2;
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
        
        % Checks if the first optional argument is DATATYPE
        elseif strcmp(varargin{4},'velocity')||strcmp(varargin{4},'ps')
            datatype=varargin{4};
            % Checks that the second, third, fourth, and fifth optional
            % arguments are METHOD, NUM_X, NUM_Y, and NUM_Z
            if strcmp(varargin{5},'interp')&&isnumeric(varargin{6})&&isnumeric(varargin{7})&&isnumeric(varargin{8})
                method=varargin{5};
                num_x=varargin{6};
                num_y=varargin{7};
                num_z=varargin{8};
                display='slice';
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
            
        % Otherwise, show help file and quit routine
        else
            help ray_plotmodel
            return;
        end
        
               
        
        % Checks if six optional arguments are given
        case 9
        r = varargin{1};
        direction=varargin{2};
        line=varargin{3};

        % Checks if the first optional argument is METHOD
        if strcmp(varargin{4},'interp')
            method=varargin{4};
            % Checks that the second, third, and fourth optional arguments
            % are NUM_X, NUM_Y, and NUM_Z
            if isnumeric(varargin{5})&&isnumeric(varargin{6})&&isnumeric(varargin{7})
                num_x=varargin{5};
                num_y=varargin{6};
                num_z=varargin{7};
                % Checks if the fifth optional argument is DISPLAY
                if strcmp(varargin{8},'slice')||strcmp(varargin{8},'contour')
                    display=varargin{8};
                    % Checks if the sixth optional argument is INTERVAL
                    if isnumeric(varargin{9})
                        interval=varargin{9};
                        datatype='velocity';
                    % Checks if the sixth optional argument is DATATYPE
                    elseif strcmp(varargin{9},'velocity')||strcmp(varargin{9},'ps')
                        datatype=varargin{9};
                        interval=0.2;
                    % Otherwise, show help file and quit routine
                    else
                        help ray_plotmodel
                        return;
                    end
                % Checks if the fifth optional argument is DATATYPE
                elseif strcmp(varargin{8},'velocity')||strcmp(varargin{8},'ps')
                    datatype=varargin{8};
                    % Checks if the sixth optional argument is DISPLAY
                    if strcmp(varargin{9},'slice')||strcmp(varargin{9},'contour')
                        display=varargin{9};
                        interval=0.2;
                    % Otherwise, show help file and quit routine
                    else
                        help ray_plotmodel
                        return;
                    end
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Otherwise, show help file and quit routine 
            else
                help ray_plotmodel
                return;
            end
        % Checks if the first optional argument is DATATYPE   
        elseif strcmp(varargin{4},'velocity')||strcmp(varargin{4},'ps')
            datatype=varargin{4};
            % Checks if the second optional argument is DISPLAY
            if strcmp(varargin{5},'slice')||strcmp(varargin{5},'contour')
                display=varargin{5};
                % Checks that third, fourth, fifth, and sixth optional
                % arguments are METHOD, NUM_X, NUM_Y, and NUM_Z
                if strcmp(varargin{6},'interp')&&isnumeric(varargin{7})&&isnumeric(varargin{8})&&isnumeric(varargin{9})
                    method=varargin{6};
                    num_x=varargin{7};
                    num_y=varargin{8};
                    num_z=varargin{9};
                    interval=0.2;
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Checks if the second, third, fourth, and fifth optional 
            % arguments are METHOD, NUM_X, NUM_Y, and NUM_Z
            elseif strcmp(varargin{5},'interp')&&isnumeric(varargin{6})&&isnumeric(varargin{7})&&isnumeric(varargin{8})
                method=varargin{5};
                num_x=varargin{6};
                num_y=varargin{7};
                num_z=varargin{8};
                % Checks if the sixth optional argument is DISPLAY
                if strcmp(varargin{9},'slice')||strcmp(varargin{9},'contour')
                    display=varargin{9};
                    interval=0.2;
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
        % Checks if first optional argument is DISPLAY    
        elseif strcmp(varargin{4},'slice')||strcmp(varargin{4},'contour')
            display=varargin{4};
            % Checks if second optional argument is INTERVAL
            if isnumeric(varargin{5})
                interval=varargin{5};
                % Checks if third, fourth, fifth, and sixth optional
                % arguments are METHOD, NUM_X, NUM_Y, and NUM_Z
                if strcmp(varargin{6},'interp')&&isnumeric(varargin{7})&&isnumeric(varargin{8})&&isnumeric(varargin{9})
                    method=varargin{6};
                    num_x=varargin{7};
                    num_y=varargin{8};
                    num_z=varargin{9};
                    datatype='velocity';
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Checks if second optional argument is DATATYPE    
            elseif strcmp(varargin{5},'velocity')||strcmp(varargin{5},'ps')
                datatype=varargin{5};
                % Checks if third, fourth, fifth, and sixth optional
                % arguments are METHOD, NUM_X, NUM_Y, and NUM_Z
                if strcmp(varargin{6},'interp')&&isnumeric(varargin{7})&&isnumeric(varargin{8})&&isnumeric(varargin{9})
                    method=varargin{6};
                    num_x=varargin{7};
                    num_y=varargin{8};
                    num_z=varargin{9};
                    interval=0.2;
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Checks if second, third, fourth and fifth optional arguments
            % are METHOD, NUM_X, NUM_Y, and NUM_Z
            elseif strcmp(varargin{5},'interp')&&isnumeric(varargin{6})&&isnumeric(varargin{7})&&isnumeric(varargin{8})
                method=varargin{5};
                num_x=varargin{6};
                num_y=varargin{7};
                num_z=varargin{8};
                % Checks if the sixth optional argument is DATATYPE
                if strcmp(varargin{9},'velocity')||strcmp(varargin{9},'ps')
                    datatype=varargin{9};
                    interval=0.2;
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
        % Otherwise, show help file and quit routine    
        else
            help ray_plotmodel
            return;
        end
        
        

        
        % Checks if all seven optional arguments are given
        case 10
        r=varargin{1};
        direction=varargin{2};
        line=varargin{3};
        
        % Checks if first, second, third, and fourth optional arguments are 
        % METHOD, NUM_X, NUM_Y, and NUM_Z
        if strcmp(varargin{4},'interp')&&isnumeric(varargin{5})&&isnumeric(varargin{6})&&isnumeric(varargin{7})
            method=varargin{4};
            num_x=varargin{5};
            num_y=varargin{6};
            num_z=varargin{7};
            % Checks if fifth and sixth optional arguments are DISPLAY and
            % INTERVAL
            if strcmp(varargin{8},'contour')&&isnumeric(varargin{9})
                display=varargin{8};
                interval=varargin{9};
                % Checks if seventh optional argument is DATATYPE
                if strcmp(varargin{10},'velocity')||strcmp(varargin{10},'ps')
                    datatype=varargin{10};
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Checks if fifth optional argument is DATATYPE    
            elseif strcmp(varargin{8},'velocity')||strcmp(varargin{8},'ps')
                datatype=varargin{8};
                % Checks if sixth and seventh optional arguments are
                % DISPLAY and INTERVAL
                if strcmp(varargin{9},'contour')&&isnumeric(varargin{10})
                    display=varargin{9};
                    interval=varargin{10};
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
            
        % Checks if first and second optional arguments are DISPLAY and INTERVAL
        elseif strcmp(varargin{4},'contour')&&isnumeric(varargin{5})
            display=varargin{4};
            interval=varargin{5};
            % Checks if third, fourth, fifth, and sixth optional arguments
            % are METHOD, NUM_X, NUM_Y, and NUM_Z
            if strcmp(varargin{6},'interp')&&isnumeric(varargin{7})&&isnumeric(varargin{8})&&isnumeric(varargin{9})
                method=varargin{6};
                num_x=varargin{7};
                num_y=varargin{8};
                num_z=varargin{9};
                % Checks if seventh optional argument is DATATYPE
                if strcmp(varargin{10},'velocity')||strcmp(varargin{10},'ps')
                    datatype=varargin{10};
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Checks if third optional argument is DATATYPE    
            elseif strcmp(varargin{6},'velocity')||strcmp(varargin{6},'ps')
                datatype=varargin{6};
                % Checks if fourth, fifth, sixth, and seventh optional arguments
                % are METHOD, NUM_X, NUM_Y, and NUM_Z
                if strcmp(varargin{7},'interp')&&isnumeric(varargin{8})&&isnumeric(varargin{9})&&isnumeric(varargin{10})
                    method=varargin{7};
                    num_x=varargin{8};
                    num_y=varargin{9};
                    num_z=varargin{10};
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Otherwise, show help file and quit routine    
            else
                help ray_plotmodel
                return;
            end
            
        % Checks if first optional argument is DATATYPE
        elseif strcmp(varargin{4},'velocity')||strcmp(varargin{4},'ps')
            datatype=varargin{4};
            % Checks if second and third optional arguments are DISPLAY and
            % INVERVAL
            if strcmp(varargin{5},'contour')&&isnumeric(varargin{6})
                display=varargin{5};
                interval=varargin{6};
                % Checks if fourth, fifth, sixth, and seventh optional
                % arguments are METHOD, NUM_X, NUM_Y, and NUM_Z
                if strcmp(varargin{7},'interp')&&isnumeric(varargin{8})&&isnumeric(varargin{9})&&isnumeric(varargin{10})
                    method=varargin{7};
                    num_x=varargin{8};
                    num_y=varargin{9};
                    num_z=varargin{10};
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Checks if second, third, fourth, and fifth optional arguments
            % are METHOD, NUM_X, NUM_Y, and NUM_Z
            elseif strcmp(varargin{5},'interp')&&isnumeric(varargin{6})&&isnumeric(varargin{7})&&isnumeric(varargin{8})
                method=varargin{5};
                num_x=varargin{6};
                num_y=varargin{7};
                num_z=varargin{8};
                % Checks if sixth and seventh optional arguments are
                % DISPLAY and INVERVAL
                if strcmp(varargin{9},'contour')&&isnumeric(varargin{10})
                    display=varargin{9};
                    interval=varargin{10};
                % Otherwise, show help file and quit routine
                else
                    help ray_plotmodel
                    return;
                end
            % Otherwise, show help file and quit routine   
            else
                help ray_plotmodel
                return;
            end
        % Otherwise, show help file and quit routine    
        else
            help ray_plotmodel
            return;
        end
            
        
        % Otherwise, show help file and quit routine
        otherwise
        help ray_plotmodel
        return;
    end
    
    if strcmp(method,'facet')
        [X,Y,Z,w]=make_griddata_facet(r,datatype);
    else
        [X,Y,Z,w]=make_griddata_interp(r,datatype,num_x,num_y,num_z);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checks if no optional arguments are given
elseif nargin==6
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    w=varargin{4};
    direction=varargin{5};
    line=varargin{6};
    method = 'interp';
    display='slice';
    
% ADD DATATYPE OPTION FROM HERE AND BELOW 

% Checks if one optional argument is given
elseif nargin==7
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    w=varargin{4};
    direction=varargin{5};
    line=varargin{6};
    % Checks if first optional argument is DISPLAY
    if (strcmp(varargin{7},'slice')||strcmp(varargin{7},'contour'))
        display=varargin{7};
        if strcmp(display,'contour')
            interval=0.2;
        end
        method='interp';
    % Checks if second optional argument is METHOD    
    elseif (strcmp(varargin{7},'interp')||strcmp(varargin{7},'facet'))
        method=varargin{7};
        display='slice';
    % Otherwise, show help file and quit routine
    else
        help ray_plotmodel
        return;
    end
    

% Checks if two optional arguments are given
elseif nargin==8
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    w=varargin{4};
    direction=varargin{5};
    line=varargin{6};
    % Checks if first optional argument is DISPLAY
    if strcmp(varargin{7},'slice')||strcmp(varargin{7},'contour')
        display=varargin{7};
        % Checks if second optional argument is INTERVAL
        if isnumeric(varargin{8})
            interval=varargin{8};
        % Checks if second optional argument is METHOD
        elseif strcmp(varargin{8},'interp')||strcmp(varargin{8},'facet')
            method=varargin{8};
            interval=0.2;
        % Otherwise, show help file and quit routine
        else
            help ray_plotmodel
            return;
        end
    % Checks if first optional argument is METHOD
    elseif strcmp(varargin{7},'interp')||strcmp(varargin{7},'facet')
        method=varargin{7};
        % Checks of second optional argument is DISPLAY
        if strcmp(varargin{8},'slice')||strcmp(varargin{8},'contour')
            display=varargin{8};
            interval=0.2;
        % Otherwise, show help file and quit routine
        else
            help ray_plotmodel
            return;
        end
    % Otherwise, show help file and quit routine
    else
        help ray_plotmodel
        return;
    end
    
% Checks if all three optional arguments are given    
elseif nargin==9
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    w=varargin{4};
    direction=varargin{5};
    line=varargin{6};
    % Checks if first optional argument is DISPLAY
    if strcmp(varargin{7},'slice')||strcmp(varargin{7},'contour')
        display=varargin{7};
        % Checks if second optional argument is INTERVAL
        if isnumeric(varargin{8})
            interval=varargin{8};
            % Checks if third optional argument is METHOD
            if strcmp(varargin{9},'interp')||strcmp(varargin{9},'facet')
                method=varargin{9};
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end        
        % Otherwise, show help file and quit routine
        else
            help ray_plotmodel
            return;
        end
    % Checks if first optional argument is METHOD
    elseif strcmp(varargin{7},'interp')||strcmp(varargin{7},'facet')
        method=varargin{7};
        % Checks if second optional argument is DISPLAY
        if strcmp(varargin{8},'slice')||strcmp(varargin{8},'contour')
            display=varargin{8};
            % Checks if third optional argument is INTERVAL
            if isnumeric(varargin{9})
                interval=varargin{9};
            % Otherwise, show help file and quit routine
            else
                help ray_plotmodel
                return;
            end
        % Otherwise, show help file and quit routine
        else
            help ray_plotmodel
            return;
        end
    % Otherwise, show help file and quit routine
    else
        help ray_plotmodel
        return;
    end
                    
        
    
% Otherwise, show help file and quit routine
else
    help ray_plotmodel
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creates a slice plot
if strcmp(display,'slice')
    hPlot=ray_plotslice(X,Y,Z,w,direction,line,method);

% Creates a contour plot
elseif strcmp(display,'contour')
    hPlot=ray_plotcontour(X,Y,Z,w,direction,line,interval);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL FUNCTION OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargout
    case 0
        return
    case 1
        varargout{1}=hPlot;
    case 4
        varargout{1}=X;
        varargout{2}=Y;
        varargout{3}=Z;
        varargout{4}=w;
    case 5
        varargout{1}=X;
        varargout{2}=Y;
        varargout{3}=Z;
        varargout{4}=w;
        varargout{5}=hPlot;
    otherwise
        disp('Error: Invalid number of output arguments');
        return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to grid the data when 'interp' option selected
function [X,Y,Z,w]=make_griddata_interp(r,datatype,num_x,num_y,num_z)
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
if strcmp(datatype,'velocity')
    w = griddata3(r.modelxyz(:,1),r.modelxyz(:,2),r.modelxyz(:,3),r.modelxyz(:,4),X,Y,Z,'linear');
elseif strcmp(datatype,'ps')
    w = griddata3(r.modelxyz(:,1),r.modelxyz(:,2),r.modelxyz(:,3),r.modelxyz(:,7),X,Y,Z,'linear');
end


% Function to grid the data when 'facet' option selected
function [X,Y,Z,w]=make_griddata_facet(r,datatype)
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

if strcmp(datatype,'velocity')
    w = griddata3(r.modelxyz(:,1),r.modelxyz(:,2),r.modelxyz(:,3),r.modelxyz(:,4),X,Y,Z,'linear');
elseif strcmp(datatype,'ps')
    w = griddata3(r.modelxyz(:,1),r.modelxyz(:,2),r.modelxyz(:,3),r.modelxyz(:,7),X,Y,Z,'linear');
end