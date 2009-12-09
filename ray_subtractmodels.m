function r = ray_subtractmodels(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_subtractmodels
%
% This function takes two RAYPACK3D structures and subtracts the starting 
% modelijk and modelxyz velocity and hits values from the final modelijk
% and modelxyz velocity and hits values node-by-node. Note that this
% function requires that the two structures be gridded the same. The output
% can be in either absolute change or percent change from the
% starting model. The resulting structure only contains the modeldims,
% modelijk, and modelxyz fields with all the other fields cleared (as they
% no longer have a logical meaning).  This function subtracts both the
% velocity field as well as the q field.
% 
% Usage: r = ray_subtractmodels(final,starting,[percent])
%
% Required Inputs:
%   final:          The RAYPACK3D structure to be subtracted from
%
%   starting:       The initial RAYPACK3D structure to subtract from the
%                   final structure
%
% Optional Inputs:
%   percent:        A boolean value (true or false) signifying if the
%                   output should be in absolute velocity (default, false),
%                   or in percent change from the starting model (true)
%
% Output:
%   r:              A RAYPACK3D structure where the r.modelijk and
%                   r.modelxyz velocity and hits values are the difference
%                   of the final model minus the starting model
%
%
% Author:
% Matt Gardine
% February 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 2
        final=varargin{1};
        starting=varargin{2};
        percent=false;
        
    case 3
        final=varargin{1};
        starting=varargin{2};
        percent=varargin{3};
        
    otherwise
        help ray_subtractmodels
        return;
end

r.hypocenter=[];
r.station=[];
r.traveltime=[];
r.residual=[];
r.rayid=[];
r.raycoord=[];
r.tetra=[];

r.modeldims=final.modeldims;
r.modelijk(:,1)=final.modelijk(:,1);
r.modelijk(:,2)=final.modelijk(:,2);
r.modelijk(:,3)=final.modelijk(:,3);
r.modelijk(:,5)=final.modelijk(:,5);
r.modelijk(:,6)=final.modelijk(:,6);

r.modelxyz(:,1)=final.modelxyz(:,1);
r.modelxyz(:,2)=final.modelxyz(:,2);
r.modelxyz(:,3)=final.modelxyz(:,3);
r.modelxyz(:,5)=final.modelxyz(:,5);
r.modelxyz(:,6)=final.modelxyz(:,6);

if percent
    for i=1:length(r.modelijk(:,4))
        r.modelijk(i,4)=(final.modelijk(i,4)-starting.modelijk(i,4))*100/starting.modelijk(i,4);
        r.modelxyz(i,4)=(final.modelxyz(i,4)-starting.modelxyz(i,4))*100/starting.modelxyz(i,4);
        r.modelijk(i,6)=(final.modelijk(i,6)-starting.modelijk(i,6));
        r.modelxyz(i,6)=(final.modelxyz(i,6)-starting.modelxyz(i,6));
        r.modelijk(i,7)=(final.modelijk(i,7)-starting.modelijk(i,7))*100/starting.modelijk(i,7);
        r.modelxyz(i,7)=(final.modelxyz(i,7)-starting.modelxyz(i,7))*100/starting.modelxyz(i,7);
    end
    
else
    for i=1:length(r.modelijk(:,4))
        r.modelijk(i,4)=(final.modelijk(i,4)-starting.modelijk(i,4));
        r.modelxyz(i,4)=(final.modelxyz(i,4)-starting.modelxyz(i,4));
        r.modelijk(i,6)=(final.modelijk(i,6)-starting.modelijk(i,6));
        r.modelxyz(i,6)=(final.modelxyz(i,6)-starting.modelxyz(i,6));
        r.modelijk(i,7)=(final.modelijk(i,7)-starting.modelijk(i,7));
        r.modelxyz(i,7)=(final.modelxyz(i,7)-starting.modelxyz(i,7));
    end
end