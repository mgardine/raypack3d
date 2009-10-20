%% RAYPACK3D Cookbook
% This is a basic introduction to using the RAYPACK3D toolbox. The purpose
% of this toolbox is to provide a Matlab interface for the raytrace3d
% package, written by William Menke (Menke 2005, 
% ftp://ftp.ldeo.columbia.edu/pub/menke/raytrace3d.tar.Z). Raytrace3d is
% command-line based code, written in C, for tracing rays through a three
% dimensional isotropic structure and performing seismic traveltime
% inversions. The inputs and outputs for this code are specifically
% formatted plain text files.  RAYPACK3D is designed to separate the user
% from having directly manipulate these files and run assorted raytrace3d
% commands from a Matlab environment.  RAYPACK3D also includes functions to
% plot inversion results as produced by raytrace3d.
%
% The purpose of this cookbook is twofold; first, to list the functions
% that are available to the user, and second, to show how some of the
% functions can be used to produce the desired results. This is not meant
% to completely document each function, that task is easily accomplished by
% typing 'help _function_' in the Matlab command window. It is also not
% meant to be a self-contained demo in the true "cookbook" fashion. Many of
% the functions require external calls to either raytrace3d or Antelope or
% both, and the results of such calls cannot be easily emulated using
% Matlab alone.


%% Functions available in RAYPACK3D
% This is a list of all of the callable functions in the RAYPACK3D toolbox
% and their purpose.
%
% *RAYPACK3D structure functions*
%
% * *ray_loaddata* : Loads the data from raytrace3d text files into a RAYPACK3D structure
%
% * *ray_subtractmodels*: Calculates the P-wave velocity differences
% between two models at each grid node
%
% *RAYPACK3D forward-modeling functions*
%
% * *ray_shoot_predicted* :   Interface to *_raytrace3d source_to_receivers_*
%
% *RAYPACK3D inversion functions*
%
% * *ray_params_gui* : A graphical interface for generating the
% model_pars_file requred by the *_raytrace3d invert_* routine
%
% * *ray_invert_db* : Interface to *_raytrace3d invert_* using arrivals
% from an Antelope database
%
% * *ray_invert_ttfile* : Interface to *_raytrace3d invert_* using arrivals
% calculated by *_raytrace3d source_to_receivers_* (or *ray_shoot_predicted*)
%
% * *ray_locate* : Interface to *_raytrace3d locate_* using arrivals given
% by an Antelope database and a raytrace3d-compatable velocity model file
%
% *RAYPACK3D plotting functions*
%
% * *ray_plotmodel* : Creates a volumetric slice plot through a velocity
% model given by a RAYPACK3D structure
%
% * *ray_plothits* : Creates a volumetric slice plot through a model given
% by a RAYPACK3D structure; plots the number of rays passing through each
% node rather than velocity at each node
%
% * *ray_plotrays* : Plots the seismic rays given by a RAYPACK3D structure
%
% *RAYPACK3D utilities*
%
% * *ray_cookbook* : The code for this cookbook
%
% * *ray_viewcookbook* : Publishes and opens the cookbook for viewing
%
% * *ray_defaults* : A user-editable file containing the default parameters
% for reference latitude, longitude, and projection type
%
% * *ray_french* : Generates a blue-to-white-to-red colormap useful for
% plotting tomographic results
%
% * *ray_lostrays* : Creates a list of origin-station pairs that exist in
% an Antelope database but cannot be connected by the *_raytrace3d_*
% routine
%
% *RAYPACK3D internal functions*
% These functions are called internally by the other RAYPACK3D functions
% described above.  They should not be modified by the user for any reason,
% as doing so may result in unintended behavior in higher-level functions.
%
% * *ray_latlon2xyz* : Converts lat/lon/elev coordinates to XYZ
% (north-east-down positive) coordinates using a spherical-earth projection
%
% * *ray_latlon2xyz_flat* : Converts lat/lon/elevation coordinates to XYZ
% (north-east-down positive) coordinates using a flat-earth projection
%
% * *ray_xyz2latlon* : Converts XYZ (north-east-down positive) coordinates
% to lat/lon/elevation coordinates using a spherical-earth projection
%
% * *ray_xyz2latlon_flat* : Converts XYZ (north-east-down positive)
% coordinates to lat/lon/elevation coordinates using a flat-earth
% projection
%
% * *ray_plotcontour* : Called by *ray_plotmodel*, creates a contour plot
%
% * *ray_plotslice* : Called by *ray_plotmodel*, creates a slice plot


%% Loading data into a RAYPACK3D structure
% Loading data is accomplished using the ray_loaddata command. It can
% currently load the output text files from ray_shoot_predicted,
% ray_locate, ray_invert_db, and ray_invert_ttfile.  The first input *MUST*
% be a model file, but the others are optional and can done in any order.
%
%>> my_ray=ray_loaddata('my_model_file.model','my_hits_file.model.hits','my_tt_file.tt');

%%
% Or,

%>> my_ray=ray_loaddata('my_model_file.model','my_locations.loc');

%%
% etc. The resulting structure has 10 fields.
%
% * *hypocenter* : Contains X,Y,Z, and origin time data for earthquakes 
% associated with the data set. Data loaded from a traveltime file (output
% from ray_shoot_predicted) or a location file (output from ray_locate). If
% both datasets exist, then the hypocenter field is loaded from the
% traveltime file.
%
% * *station* : Contains the X,Y,Z, and arrival time data for each arrival
% associated with the data set. Data is loaded from the traveltime file
% (output from ray_shoot_predicted).
%
% * *traveltime* : Arrival time - origin time, given in seconds. This field
% is derived from the traveltime file (output from ray_shoot_predicted).
%
% * *residual* : Observed arrival time - predicted arrival time, given in
% seconds. This field is derived from the traveltime file (output from
% ray_shoot_predicted).
%
% * *rayid* : A unique id number for each ray connecting an origin to a
% station. This field is loaded from the rayfile (output from
% ray_shoot_predicted).
%
% * *raycoord* : A series of X,Y,Z coordinates for each ray connecting an
% origin to a station. This field is loaded from the rayfile (output from
% ray_shoot_predicted).
%
% * *tetra* : An index number corresponding to which tetrahedron the ray is
% passing through at each raycoord point. This field is loaded from the
% rayfile (output from ray_shoot_predicted).
%
% * *modeldims* : The number of X, Y, and Z nodes in the model. This field
% is loaded from the required model file given as the first input.
%
% * *modelijk* : The given model shown in index notation (i,j,k), the
% velocity at each node, the interface number of the node, and the number
% of rays that pass through the node (loaded from model.hits file, output
% from ray_invert_db and ray_invert_ttfile). If no model.hits file is
% given, then the final column is set to zero.
%
% * *modelxyz* : The given model shown in x,y,z coordinates (km), the
% velocity at each node, the interface number of the node, and the number
% of rays that pass through the node (loaded from model.hits file, output
% from ray_invert_db and ray_invert_ttfile). If no model.hits file is
% given, then the final column is set to zero.

%% Subtracting models
% Often it is desireable to subtract two models from each other in order to
% view results as pertubations rather than absolute values. RAYPACK3D
% accomplishes this with the *ray_subtractmodels* function. The function
% takes two RAYPACK3D structures, a "final" structure and a "starting"
% structure, and subtracts the starting model from the final model 
% node-by-node. To do this properly, the models need to have the exact same
% geometry, although *ray_subtractmodels* does not explicitly check for it.
% By default, the function returns a new structure with absolute velocity
% change and the difference in the number of rays passing through each
% node.

%>>result = ray_subtractmodels('./final_structure','./starting_structure')

%%
% This function can also return the results in percent change from starting
% model by adding a third input 'true':

%>>results = ray_subtractmodels('./final_structure','./starting_structure','true')

%%
% One important note is that all of the fields in the RAYPACK3D structure
% output from this function are cleared except for modeldims, modelijk,
% and modelxyz.

%% Forward-modeling with RAYPACK3D
% Given an Antelope database with origin and station information, as well
% as a raytrace3d model file, RAYPACK3D has the ability to generate
% predicted traveltimes and rays that connect origin-station pairs. This is
% done by the *ray_shoot_predicted* function. This function requires an
% Antelope database with site, origin, event, assoc, and arrival tables. An
% optional third argument allows for subsetting the database by inputing a
% valid Datascope argument. To calculate predicted traveltimes using
% origins in a database _my_database_ and a model file _my_model_:

%>>ray_shoot_predicted('./my_database','./my_model');

%% 
% The output would be two files, a predicted traveltime file called
% _pred_tt.tt_ and a ray file called _ray_file.rays_. Alternatively, say we
% wanted to use the same database and model file, but only wanted arrivals
% at stations that are within 100 km from the origin.

%>>ray_shoot_predicted('./my_database','./my_model','delta<(100/111.1)');

%% Making a model_pars file
% One required input for all ray_invert routines is a model_pars file,
% which contains information about each node in the model file. By properly
% using this file, the user can specify which nodes much change together
% and how much weight to give measurements at each node. This file can also
% give some control over smoothing in the model. While this is a fairly
% straight-forward text file, RAYPACK3D provides a GUI for creating and
% manipulating this file graphically through the *ray_params_gui* function.
% This function requires a RAYPACK3D structure as input; it will ask for a
% model file and a model_pars file (if one already exists for editing,
% otherwise it will create a new file). The user can then click on the
% nodes in each interface to create the desired groupings.

%% P-wave tomography with RAYPACK3D
% RAYPACK3D has the ability to do tomographic inversions using P-wave
% arrival data from two sources. The first source is via a traveltime file
% such as one that is output from *ray_shoot_predicted*. This is done
% through the *ray_invert_ttfile* routine. 

%>>ray_invert_ttfile('./my_tt_file','./my_model','./my_model_pars',50)

%%
% This command would take the files _my_tt_file_, _my_model_, and
% _my_model_pars_, and run the inversion using a damping value of 50. 
%
% The second source is via an Antelope database containing an origin, event, 
% site, affiliation, arrival, and assoc table. This is done through the
% *ray_invert_db* routine. 

%>>ray_invert_db('./my_database','./my_model','./my_model_pars',50)

%%
% Again, this command would take the arrival and origin information from the
% _my_database_ database and run the inversion through _my_model_, with
% _my_model_pars_ file and a damping value of 50. Additionally,
% *ray_invert_db* allows for two optional inputs. The first is a Datascope
% argument for subsetting the database, such as:

%>>ray_invert_db('./my_database','./my_model','./my_model_pars',50,'net=~/ZA/')

%%
% would run the same inversion as before but subset the database to only
% include arrivals appearing on stations with a network code of ZA. The
% second optional input allows for origin information to be taken from the 
% output of *ray_locate* rather than from the origin table. This option is 
% especially useful in an iterative inversion process where origins are 
% relocated using a model generated by a previous run of *ray_invert_db*.

%>>ray_invert_db('./my_database','./my_model','./my_model_pars',50,'./my_locations')

%%
% In both *ray_invert_ttfile* and *ray_invert_db*, the output is a new model
% file called _inversion_model.model_ and a new model hits file (containing
% the number of rays passing through each grid node rather than velocities
% at each node) called _inversion_model.model.hits_. These files are
% easily loaded into a RAYPACK3D structure via the *ray_loaddata* command.

%% Locating earthquakes with RAYPACK3D
% RAYPACK3D provides an interface to the simple earthquake location routine
% built into raytrace3d through the *ray_location* function. This function
% calls _raytrace3d_ _locate_ and uses arrival information from an Antelope
% database. This routine uses Geiger's method for locating earthquakes and
% assumes a standard P-to-S velocity ratio of 1.76. Inputs to the function
% are the path to an Antelope database _my_database_, a model file
% _my_model_, a tolerance value, in seconds, below which the iterations
% will stop (otherwise the program will run 5 iterations, set to 0.5 seconds
% in the example), and a damping value (set to 1 in the example). An
% optional input is a Datascope argument for subsetting the database (not
% shown).

%>>ray_locate('./my_database','./my_model',0.5,1)

%% Plotting with RAYPACK3D
% Viewing volumetric results is often a complicated process. RAYPACK3D
% provides a basic framework for plotting Raytrace3d results using the
% RAYAPCK3D structure as an input.
%
% *Plotting rays*
%
% Rays connecting source to receiver can be plotted using the
% *ray_plotrays* command. This function provides two options - the user can
% simply give a RAYPACK3D structure as input to plot all of the rays loaded
% into the structure:

%>>ray_plotrays('./my_structure')

%%
% Alternatively, a latitude and longitude coordinate can be given
% corresponding to the coordinates of a particular station, and
% *ray_plotrays* will only plot the rays from all of the sources to the
% single receiver:

%>>ray_plotrays('./my_structure',19.65,-104.25)

%%
% where 19.65 is the latitude of a station and -104.25 is the longitude.
%
% *Plotting Volumetric Slices*
%
% RAYPACK3D's *ray_plotmodel* function is a highly customizable command
% that creates both volumetric slice plots as well as contour plots using a
% RAYPACK3D structure as the first input. The second input is used to
% determine what slice plane to plot. The options are 'xy' for a depth
% slice in the x-y plane, 'xz' for a cross-section in the x-z plane, and 
% 'yz' for a cross-section in the y-z plane. The third required input is a
% numerical value at which to fix the plane at in the third axis. For
% example, to create a depth slice at a depth of 20 km using a RAYPACK3D 
% structure called 'my_structure':

%>>ray_plotmodel('./my_structure','xy',20)

%%
% However, this command has many additional options. The first is the
% method for plotting. Options are 'interp' to re-sample the grid into
% equally spaced nodes and interpolate the values (necessary for grids with
% dipping layers), and 'facet' which plots the nodes as they are. By
% default, ray_plotmodel uses the 'interp' option with a 50x50x50 grid. To
% change the number of nodes to 100 x-nodes, 50 y-nodes, and 20 z-nodes:

%>>ray_plotmodel('./my_structure','xy',20,'interp',100,50,20)

%%
% Or, to plot using the 'facet' option instead of 'interp':

%>>ray_plotmodel('./my_structure','xy',20,'facet')

%%
% In the process of making a slice plot, ray_plotmodel calls the Matlab
% function griddata3 in order to make a three-dimensional array
% representing the model. This is by far the most time-consuming process in
% ray_plotmodel, and so minimizing the number of times that the arrays are
% created is often desirable. To do this, ray_plotmodel has an option to
% output the four arrays created by griddata3:

%>>[X,Y,Z,w]=ray_plotmodel('./my_structure','xy',20,'facet')

%%
% These output arrays can then be used as an input to ray_plotmodel
% (instead of a RAYPACK3D structure), in order to plot another slice
% through the model.

%>>ray_plotmodel(X,Y,Z,w,'yz',100,'facet')

%%
% would plot a cross-section in the y-z plane with x fixed at 100 km, using
% the same data as in the prior command, but it would plot much faster
% because the function call to griddata3 is not necessary.
%
% A final optional output is the Matlab graphics handle for the given plot:

%>>[my_handle]=ray_plotmodel('./my_structure','xy',20,'facet')

%%
% Or, combined with the array output:

%>>[X,Y,Z,w,my_handle]=ray_plotmodel('./my_structure','xy',20,'facet')

%%
% The graphics handle is useful for low-level manipulation of Matlab
% figures, details for which will not be given here but can be found in the
% Matlab documentation.
%
% *Plotting Contour Slices*
%
% Instead of shaded slices, ray_plotmodel can also plot contour slices. The
% syntax is almost the same, except the 'display' option should be
% 'contour' instead of 'slice':

%>>ray_plotmodel('./my_structure','xy',20,'contour')

%%
% All of the options shown for making slice plots also work for contour
% plots. In addition, an extra option is the interval value, which is used
% to specify what the contour interval should be for the plot. For example,
% to make a plot with contour lines every 0.5 km/s (as opposed to the
% default value of 0.2 km/s):

%>>ray_plotmodel('./my_structure','xy',20,'contour',0.5)

%% Advanced Scripting Options with RAYPACK3D
% Many of the functions in RAYPACK3D can be combined to form more advanced
% routines. A common task is to invert for the velocity model, relocate the
% origins using the new velocity model, and invert the data again. This is
% easily accomplished with a script such as:

%>>ray_locate('/home/mgardine/testing/NEW/slab_events/used_events','./1d_starting.model',0.5,1);
%>>system('mv location_output.tsv location_orig.tsv');
%>>ray_invert_db('/home/mgardine/testing/NEW/slab_events/used_events','./1d_starting.model','./invert_mps.txt',50,'./location_orig.tsv')
%>>system('mv inversion_model.model inversion1.model');
%>>system('mv inversion_model.model.hits inversion1.model.hits');
%>>system('mv invert_output.log invert1.log');
%>>ray_locate('/home/mgardine/testing/NEW/slab_events/used_events','inversion1.model',0.5,1);
%>>system('mv location_output.tsv location1.tsv');
%
%>>for i=1:9
%    ray_invert_db('/home/mgardine/testing/NEW/slab_events/used_events',['./inversion' num2str(i) '.model'],...
%        './invert_mps.txt',50,['./location' num2str(i) '.tsv']);
%    runme1=['mv inversion_model.model inversion' num2str(i+1) '.model'];
%    runme2=['mv inversion_model.model.hits inversion' num2str(i+1) '.model.hits'];
%    runme3=['mv invert_output.log invert' num2str(i+1) '.log'];
%    system(runme1);
%    system(runme2);
%    system(runme3);
%    
%    ray_locate('/home/mgardine/testing/NEW/slab_events/used_events',['./inversion' num2str(i+1) '.model'],0.5,1);
%    runme4=['mv location_output.tsv location' num2str(i+1) '.tsv'];
%    system(runme4);
% end

%%
% Another common script is one that combines multiple plotting functions to
% produce one figure. The following script takes a RAYPACK3D structure
% called 'my_structure' and plot a volumetric slice at 25km depth and then
% overlays it with a contour plot with simple black lines of the same data 
% to enhance the clarity of the figure.

%>>[X,Y,Z,w]=ray_plotmodel(my_structure,'xy',25,'facet');
%>>shading interp
%>>hold on
%>>hPlot=ray_plotmodel(X,Y,Z,w,'xy',25,'contour',0.1,'facet');
%>>set(hPlot,'EdgeColor',[0 0 0]);
%>>axis equal
%>>xlim([230 320])
%>>ylim([205 295])
%>>zlim([0 65])
%>>set(gca,'FontName','Helvetica');
%>>hTitle = title('Depth 25 km');
%>>set(hTitle,'FontName','Helvetica', 'FontSize', 11, 'FontWeight','bold','units','inches');
