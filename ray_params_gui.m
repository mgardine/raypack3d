function ray_params_gui(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_params_gui
%
% This function is a graphical front-end for creating the model_pars_file 
% parameter file necessary for the "raytrace3d invert" routine, as well as 
% the ray_invert_ttfile and ray_invert_db matlab routines.  It can create a
% new file, or load a previously saved model_pars_file for editing.  With
% this GUI, the user can create new groups of grid nodes (Create button),
% add nodes to existing groups (Associate button), remove nodes from
% existing groups (Unassociate button), find out which group a node is
% associated with (Which button), copy the group structure of the current
% interface to other interfaces (Copy button), and write out the final
% model_pars_file (Write button).
%
% Usage: ray_params_gui(r)
%
% Required Inputs:
%   r:          A matlab raytrace3d structure
%
% Output (to file system):
%   A user-specified filename for the model_pars_file when the "write" 
%   button is clicked.
%
%
%
% Author:
% Matt Gardine
% March 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch nargin
    case 1
        if isstruct(varargin{1})
            r = varargin{1};
        
            [r,g,colors] = editparams(r);
            DOTRI=1;
            doint=0;
            plotint(r,g,DOTRI,doint,colors);
        else
            disp('Error: Input argument must be a raytrace structure')
        end
    otherwise
        disp('Usage: ray_params_gui(ray_struct)')
end




function [r,g,colors] = editparams(r)
 
% READ PARAMETER FILE
disp(' ');
disp('Name of inversion parameter file,');
INP=input('  (<ENTER> to create new file) : ','s');
if isempty(INP)
   g=zeros(size(r.modelijk,1),1);
   
   [tmp index]=sortrows([r.modelijk(:,1) r.modelijk(:,2) r.modelijk(:,3)],[3 2 1]);
   r.modelijk(:,:)=r.modelijk(index,:);
   r.modelxyz(:,:)=r.modelxyz(index,:);
   
else
   [filelength]=textread(INP,'%s','whitespace','\n');
   filelength=length(filelength);
   nodei=[];
   nodej=[];
   nodek=[];
   nodegroup=[];
   group=1;
   toskip=0;
   while toskip<filelength
      [howmany]=textread(INP,'%s',1,'headerlines',toskip);
      howmany=str2num(char(howmany));
      toskip=toskip+1;
      for toskip=toskip:toskip+howmany-1
         [tmpi tmpj tmpk]=textread(INP,'%d %d %d',1,'headerlines',toskip);
         nodei=[nodei ; tmpi];
         nodej=[nodej ; tmpj];
         nodek=[nodek ; tmpk];
         nodegroup=[nodegroup ; group];
      end
      toskip=toskip+1;
      group=group+1;
   end
   
      % SORT MODEL FILE INTO SPECIFIED FORMAT
   [tmp index]=sortrows([r.modelijk(:,1) r.modelijk(:,2) r.modelijk(:,3)],[3 2 1]);
   
   r.modelijk(:,:)=r.modelijk(index,:);
   r.modelxyz(:,:)=r.modelxyz(index,:);
   
   % CREATE GROUP NUMBER VECTOR
   g=zeros(size(r.modelijk,1),1);
   for nindex=1:length(nodegroup)  %INDEX THROUGH NODES
     mindex=1+nodei(nindex)+nodej(nindex)*r.modeldims(1)+nodek(nindex)*r.modeldims(1)*r.modeldims(2);
     g(mindex)=nodegroup(nindex);
   end
   clear nindex mindex
end

% CREATE PLOT AND RANDOM GROUP COLOR SCHEME (UPTO 1000 GROUPS)
figure('Position',[30 30 700 700],'Color','w') 
colors=0.5*rand(size(r.modelijk,1)+10,3)+.5;


function plotint(r,g,DOTRI,doint,colors)
% This script is run within editparam.m
% it plots a selected interface

%global r g DOTRI doint index G I J K V INT X Y Z Hint

numint=max(r.modelijk(:,5))+1;

% GET INTERFACE POINTS / SORT BY GROUP
index=[doint*r.modeldims(1)*r.modeldims(2):(doint+1)*r.modeldims(1)*r.modeldims(2)-1]'+1;
G=g(index);                    % SELECT POINTS
[junk,index]=sort(G);          % SORT BY GROUP
index=index+doint*r.modeldims(1)*r.modeldims(2);   % GRAB SORTED NODES
G=g(index);

X=r.modelxyz(index,1);
Y=r.modelxyz(index,2);


% PLOT EACH GROUP ON THIS INTERFACE
hold off;
clf;
set(gcf,'DefaultAxesFontSize',14);
set(gcf,'DefaultAxesDataAspectRatio',[1 1 1]);
set(gcf,'DefaultLineLineWidth',2);
set(gcf,'DefaultLineMarkerSize',7);
set(gcf,'DefaultAxesXDir','reverse');
hold on;
box on;
axis([min(r.modelxyz(:,1))-5 max(r.modelxyz(:,1))+5 min(r.modelxyz(:,2))-5 max(r.modelxyz(:,2))+5]);
%axis([min(r.modelxyz(:,2))-5 max(r.modelxyz(:,2))+5 min(r.modelxyz(:,1))-5 max(r.modelxyz(:,1))+5]);
 
% PLOT CONNECTING TRIANGLES
for group=1:max(G)
   tmpX=[];
   tmpY=[];
   found=find(~(G-group));
   tmpX=X(found)+0.0001*rand(size(found));
   tmpY=Y(found)+0.0001*rand(size(found));
   if length(found)>2
     if DOTRI
        tri=delaunay(tmpX,tmpY);
        for trinum=1:length(tri(:,1))
           fill(tmpX(tri(trinum,:)),tmpY(tri(trinum,:)),colors(group,:));
           plot(tmpX(tri(trinum,[1 2 3 1])),tmpY(tri(trinum,[1 2 3 1])),'ko-','MarkerFaceColor',colors(group,:));
        end 
     else
        plot(tmpX,tmpY,'ko','MarkerFaceColor',colors(group,:));
     end;
      
      
   elseif (length(found)==2) || (length(found)==1)
      plot(tmpX,tmpY,'ko-','MarkerFaceColor',colors(group,:));
   %else
   %  if length(found)==0
   %      disp(['There are no longer any points in group ' num2str(group)]);
   %  end
   end
end
plot(X,Y,'ko','LineWidth',1,'MarkerFaceColor','none');      % PLOTS NODES
xlabel('X');
ylabel('Y');
title(['Current node groups on interface ' num2str(doint)],'FontSize',16);
view(90,90);


% ADD BUTTONS
intstring=[];
for index=1:numint
    intstring=[intstring ; {['layer ' num2str(index-1)]}];
end

uicontrol('Style','popup','ListBoxTop',doint,'BackGroundColor','w','String',...
   intstring,'Position',[20 10 70 20],...
   'Callback',{@selectlayer,r,g,DOTRI,colors});
uicontrol('Style','pushbutton','BackGroundColor','w','String',...
   'unassociate','Position',[95 10 80 20],...
   'Callback',{@unassociate,r,g,DOTRI,doint,colors});
uicontrol('Style','pushbutton','BackGroundColor','w','String',...
   'associate','Position',[180 10 70 20],...
   'Callback',{@associate,r,g,DOTRI,doint,colors});
uicontrol('Style','pushbutton','BackGroundColor','w','String',...
   'create','Position',[260 10 70 20],...
   'Callback',{@create,r,g,DOTRI,doint,colors});
uicontrol('Style','pushbutton','BackGroundColor','w','String',...
   'which','Position',[340 10 70 20],...
   'Callback',{@whichgroup,r,g,doint});
uicontrol('Style','pushbutton','BackGroundColor','w','String',...
   'copy','Position',[420 10 70 20],...
   'Callback',{@copylayer,r,g,DOTRI,doint,colors});
uicontrol('Style','pushbutton','BackGroundColor','w','String',...
   'write','Position',[500 10 70 20],...
   'Callback',{@writeout,r,g});
if DOTRI==0
  DOTRILABEL='triangulate';
else
  DOTRILABEL='no triangles';
end;
uicontrol('Style','pushbutton','BackGroundColor','w','String',...
   DOTRILABEL,'Position',[580 10 70 20],...
   'Callback',{@dotri,r,g,DOTRI,doint,colors});



function selectlayer(hObject,eventdata,r,g,DOTRI,colors)

% Callback function for Hint
% Select a layer and plot it

% VAL IN LAYER NUMBER
val=get(hObject,'Value')-1;
doint=val;
plotint(r,g,DOTRI,doint,colors);

function unassociate(hObject,eventdata,r,g,DOTRI,doint,colors)

% Callback from unassociate button
% unassociate node(s) from a group


% SELECT POINT(S) 
disp(' ');
disp('Click on point(s) to remove from group')
INP=ginput;


% FIND NODE AND REMOVE
index=[doint*r.modeldims(1)*r.modeldims(2):(doint+1)*r.modeldims(1)*r.modeldims(2)-1]'+1;
tri=delaunay(r.modelxyz(index,1),r.modelxyz(index,2),{'QJ'});
for INPnum=1:length(INP(:,1))
   NODE=dsearch(r.modelxyz(:,1),r.modelxyz(:,2),tri,INP(INPnum,1),INP(INPnum,2));
   NODE=NODE+doint*r.modeldims(1)*r.modeldims(1);
   %disp(['removing NODE ' num2str(NODE)]);
   g(NODE)=0;
end



% REPLOT FIGURE
plotint(r,g,DOTRI,doint,colors);

function associate(hObject,eventdata,r,g,DOTRI,doint,colors)
% Callback from associate button
% associate node(s) with a group


% SELECT POINT(S)
disp(' ');
disp('Click on point(s) to include in group, then <return>')
INP=ginput;
disp('Select group to associate with')
INP2=ginput(1);
index=[doint*r.modeldims(1)*r.modeldims(2):(doint+1)*r.modeldims(1)*r.modeldims(2)-1]'+1;
tri=delaunay(r.modelxyz(index,1),r.modelxyz(index,2),{'QJ'});


% FIND NUMBER OF GROUP TO ASSOCIATE
NODE=dsearch(r.modelxyz(:,1),r.modelxyz(:,2),tri,INP2(1),INP2(2));
NODE=NODE+doint*r.modeldims(1)*r.modeldims(2);
newgroup=g(NODE);


% FIND NODES AND ASSOCIATE
for INPnum=1:length(INP(:,1))
   NODE=dsearch(r.modelxyz(:,1),r.modelxyz(:,2),tri,INP(INPnum,1),INP(INPnum,2));
   NODE=NODE+doint*r.modeldims(1)*r.modeldims(2);
   %disp(['Adding NODE ' num2str(NODE)]);
   g(NODE)=newgroup;
end


% REPLOT FIGURE
plotint(r,g,DOTRI,doint,colors);

function create(hObject,eventdata,r,g,DOTRI,doint,colors)

% Callback from create button
% create new group


% SELECT POINT(S) 
disp(' ');
disp('Click on point(s) to add to new group')
INP=ginput;


% FIND NODES AND CREATE GROUP
index=[doint*r.modeldims(1)*r.modeldims(2):(doint+1)*r.modeldims(1)*r.modeldims(2)-1]'+1;
tri=delaunay(r.modelxyz(index,1),r.modelxyz(index,2),{'QJ'});
newgroup=max(g)+1;
%disp(['creating new group.  Temporary label = ' num2str(newgroup)]);
for INPnum=1:length(INP(:,1))
   NODE=dsearch(r.modelxyz(:,1),r.modelxyz(:,2),tri,INP(INPnum,1),INP(INPnum,2));
   NODE=NODE+doint*r.modeldims(1)*r.modeldims(2);
   g(NODE)=newgroup;
end


% REPLOT FIGURE
plotint(r,g,DOTRI,doint,colors);

function whichgroup(hObject,eventdata,r,g,doint)
% Callback from create button
%list temporary group name


% SELECT POINT(S) 
disp(' ');
disp('Select a node')
INP=ginput(1);


% FIND group name
index=[doint*r.modeldims(1)*r.modeldims(2):(doint+1)*r.modeldims(1)*r.modeldims(2)-1]'+1;
tri=delaunay(r.modelxyz(index,1),r.modelxyz(index,2),{'QJ'});
NODE=dsearch(r.modelxyz(:,1),r.modelxyz(:,2),tri,INP(1),INP(2));
NODE=NODE+doint*r.modeldims(1)*r.modeldims(2);
if g(NODE)==0
   disp('This node is NOT in a group');  
else 
   disp('Note: group numbers are arbitrarily assigned with each execution.');
   disp(['This node is in group ' num2str(g(NODE))]);
end

function copylayer(hObject,eventdata,r,g,DOTRI,doint,colors)

% Callback from copy button of editparam.m
% Copies current layer to other layers as specified


% GET LAYERS  NUMBERS TO COPY TO
disp('This procedure copies the group layout of the current layer onto other');
disp('selected layers.  The grouping of the selected layers is overridden and');
disp('replaced.  Layers which are not selected are not affected.  Since this');
disp('procedure operates layer by layer, it may cause problems with parameter');
disp('files that includes multiple layers in one group.  Current model');
disp(['has layers 0 through' num2str(max(r.modelijk(:,5))) '.']);
disp(' ');
disp('Enter a list of space separated layer numbers (not including the current');
disp('layer) to fill with the same grouping. ex: 1 4 8 11   (<ENTER> to abort)')
INP=input(': ','s');
if isempty(INP)
  return;
end;
copylayers=str2num(INP);




% GET UNSORTED G VALUES. **CAUTION**: OTHER VARIABLES NOT REORDERED
index=[doint*r.modeldims(1)*r.modeldims(2):(doint+1)*r.modeldims(1)*r.modeldims(2)-1]'+1;
G=g(index);     


% CREATE Gshift WHICH IS GROUP LIST (G) RENUMBERED AS 1,2,3,...
Glist=[];
Gshift=1;
for Gindex=1:max(G)
  foundG=find(~(G-Gindex));   % LOCATE ALL NODES IN A GROUP
  if ~isempty(foundG)
    Glist(Gindex)=Gshift;
    Gshift=Gshift+1;
  end;
end;
Gshift=[];
for index=1:length(G)     
  if G(index)==0
     Gshift(index)=0;                  % ASSIGN ZERO IF NOT IN GROUP
  else
     Gshift(index)=Glist(G(index));    % ASSIGN RENUMBERED GROUP
  end;
end;

  
% ASSIGN G VALUES TO DIFFERENT LAYERS BASED ON Gshift
for layer=copylayers		% STEP THROUGH SELECTED LAYERS
   gmax=max(g);
   index=[layer*r.modeldims(1)*r.modeldims(2):(layer+1)*r.modeldims(1)*r.modeldims(2)-1]'+1;
   for ii=1:length(index)
      if Gshift(ii)==0
	 g(ii+min(index)-1)=0;
      else
         g(ii+min(index)-1)=Gshift(ii)+gmax;
      
      end;
   end;
end;

% REPLOT FIGURE
plotint(r,g,DOTRI,doint,colors);


function writeout(hObject,eventdata,r,g)

% Callback from unassociate button
% unassociate node(s) from a group
%

% QUERY USER
savename=input('Enter the name of the file to be saved: ','s');
INP=input(['Save file as ' savename '? (y/n)'],'s');
if (INP~='y') && (INP~='Y')
  disp('File NOT saved.')
  return
end
%savename='tmp_parameter_file';
  
% WRITE OUT FILE
fid=fopen(savename,'w');
for group=1:max(g)
   found=find(g==group);
   if ~isempty(found)
      stringout=num2str(length(found));
      fprintf(fid,'%s\n',stringout);
      for index=[found']
         stringout=[num2str(r.modelijk(index,1)) ' ' num2str(r.modelijk(index,2)) ' ' num2str(r.modelijk(index,3))];
         fprintf(fid,'%s\n',stringout);
      end
   end
end 
fclose(fid);
disp(' ');
disp(['File saved in current directory as ' savename]);

function dotri(hObject,eventdata,r,g,DOTRI,doint,colors)

% Callback function for Hdotri
% Turns triangle plotting on and off


% SET DOTRI VALUE

if DOTRI==1
  DOTRI=0;
else
  DOTRI=1;
end;

plotint(r,g,DOTRI,doint,colors);