function [lat,lon,elev]=ray_xyz2latlon_flat(x,y,z,ref_lat,ref_lon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_xyz2latlon_flat
%
% This function takes NED (North [positve x], East [positive y], Down 
% [positive z]) cartesian coordinate system and converts them to latitude & 
% longitude coordinates, assuming a flat earth. In this tranformation, 
% x, y, and z are in kilometers. For a spherical earth
% projection, use ray_xyz2latlon.
%
% Usage: [lat,lon,elev]=ray_latlon2xyz_flat(x,y,z,ref_lat,ref_lon)
%
% Inputs:
%   x:          X coordinate of data point
%   y:          Y coordinate of data point
%   z:          Z coordinate of data point
%   ref_lat:    Reference latitude where x will be equal to 0
%   ref_lon:    Reference longitude where y will be equal to 0
%
% Outputs:
%   lat:        Output latitude coordinate
%   lon:        Output longitude coordinate
%   elev:       Output elevation (in kilometers)
%
% Author:
% Matt Gardine
% February 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 6367;
ned = [x; y; z];
rotation = [0 1 0; 1 0 0; 0 0 -1];
  
enu = rotation*ned;
  
x = enu(1);
y = enu(2);
z = enu(3);

lat = ref_lat + y*180/(pi*R);
lon = ref_lon + x*180/(pi*R*cos(lat*pi/180));
elev = z;
