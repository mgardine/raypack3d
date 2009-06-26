function [ref_lat,ref_lon,projection] = ray_defaults()

% This is the default parameter file for the raytools structure
% Valid options are:
%
% ref_lat : Reference latitude for which x = 0;
%
% ref_lon : Reference longitude for which y = 0;
%
% projection : used for converting lat/lon coordinates into x,y,z coords
%   'flat' : uses a flat-earth projection
%   'spherical : uses a spherical-earth projection
%
%

ref_lat = 17.01;
ref_lon = -105.99;
projection = 'flat';