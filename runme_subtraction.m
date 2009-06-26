%load ./starting.mat
%load ./test2/test2.mat
%load ./test1/test1.mat
%load ./test3/test3.mat
%load ./test4/test4.mat
%load ./test5/test5.mat
%load ./test6/test6.mat
% load ../MODELING/test1/test1_run2.mat
% 
% 
addpath('/home/mgardine/matlab/codex/raytrace/raytrace3d/structure');

test = ray_subtractmodels(real15,original);

%junk1=ray_subtractmodels(original14,original13);
%junk2=ray_subtractmodels(real14,real13);

%test=ray_subtractmodels(junk2,junk1);

%load ./test4/test4.mat
%test2=test4;

%load checkerboard
%test=r;

figure('Color','white','Position',[20 20 1024 768])
set(gcf,'DefaultAxesFontSize',10)

subplot(3,2,1)
%ray_plotrays(run1a,19.5240,-103.6075,'b');
%hold on
[X,Y,Z,w]=ray_plotmodel(test,'xz',254,'facet');
%[X,Y,Z,w]=ray_plotmodel(test,'xy',0,'facet');
%[X,Y,Z,w]=ray_plotmodel(test,'yz',265,'facet');
xlim([244 305])
ylim([212 292])
%xlim([240 310])
%ylim([200 300])
zlim([0 60])
shading interp
caxis([-4 4])
%caxis([-1.5 0.5])
set(gca,'Position',get(gca,'OuterPosition') - 1.3 * get(gca,'TightInset')...
        * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    
    % Uses prettier fonts than the defaults
set(gca,'FontName','Helvetica');
%hTitle = title('Depth 0 km');
hTitle = title('Cross-section at Y=254');
set(hTitle,'FontName','Helvetica', 'FontSize', 11, 'FontWeight','bold','units','inches');
    % Moves the titles down slightly to better use space
P=get(hTitle,'Position');
P(2)=P(2)-.05;
set(hTitle, 'Position', P)
hColor=findobj(gcf,'tag','Colorbar');
C=get(hColor,'Position');
C(1)=C(1)-.01;
set(hColor,'Position',C);


subplot(3,2,2)
%ray_plotrays(run1a,19.5240,-103.6075,'b');
hold on
ray_plotmodel(X,Y,Z,w,'yz',265,'facet');
%ray_plotmodel(X,Y,Z,w,'xy',5,'facet');
xlim([244 305])
ylim([212 292])
%xlim([240 310])
%ylim([200 300])
zlim([0 60])
shading interp
caxis([-4 4])
%caxis([-1.5 0.5])
set(gca,'Position',get(gca,'OuterPosition') - 1.3 * get(gca,'TightInset')...
        * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    
    % Uses prettier fonts than the defaults
set(gca,'FontName','Helvetica');
hTitle = title('Cross-section at X=265');
%hTitle = title('Depth 5 km');
set(hTitle,'FontName','Helvetica', 'FontSize', 11, 'FontWeight','bold','units','inches');
P=get(hTitle,'Position');
P(2)=P(2)-.05;
set(hTitle, 'Position', P)
hColor=findobj(gcf,'tag','Colorbar');
C=get(hColor,'Position');
C{1}(1)=C{1}(1)-.01;
set(hColor(1),'Position',C{1});

subplot(3,2,3)
ray_plotmodel(X,Y,Z,w,'xy',10,'facet');
xlim([244 305])
ylim([212 292])
%xlim([240 310])
%ylim([200 300])
zlim([0 60])
shading interp
caxis([-4 4])
%caxis([-1.5 0.5])
set(gca,'Position',get(gca,'OuterPosition') - 1.3 * get(gca,'TightInset')...
        * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    
    % Uses prettier fonts than the defaults
set(gca,'FontName','Helvetica');
hTitle = title('Depth 10 km');
set(hTitle,'FontName','Helvetica', 'FontSize', 11, 'FontWeight','bold','units','inches');
P=get(hTitle,'Position');
P(2)=P(2)-.05;
set(hTitle, 'Position', P)
hColor=findobj(gcf,'tag','Colorbar');
C=get(hColor,'Position');
C{1}(1)=C{1}(1)-.01;
set(hColor(1),'Position',C{1});


subplot(3,2,4)
ray_plotmodel(X,Y,Z,w,'xy',20,'facet');
xlim([244 305])
ylim([212 292])
%xlim([240 310])
%ylim([200 300])
zlim([0 60])
shading interp
caxis([-4 4])
%caxis([-1.5 0.5])
set(gca,'Position',get(gca,'OuterPosition') - 1.3 * get(gca,'TightInset')...
        * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    
    % Uses prettier fonts than the defaults
set(gca,'FontName','Helvetica');
hTitle = title('Depth 20 km');
set(hTitle,'FontName','Helvetica', 'FontSize', 11, 'FontWeight','bold','units','inches');
P=get(hTitle,'Position');
P(2)=P(2)-.05;
set(hTitle, 'Position', P)
hColor=findobj(gcf,'tag','Colorbar');
C=get(hColor,'Position');
C{1}(1)=C{1}(1)-.01;
set(hColor(1),'Position',C{1});


subplot(3,2,5)
ray_plotmodel(X,Y,Z,w,'xy',30,'facet');
xlim([244 305])
ylim([212 292])
%xlim([240 310])
%ylim([200 300])
zlim([0 60])
shading interp
caxis([-4 4])
%caxis([-1.5 0.5])
set(gca,'Position',get(gca,'OuterPosition') - 1.3 * get(gca,'TightInset')...
        * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    
    % Uses prettier fonts than the defaults
set(gca,'FontName','Helvetica');
hTitle = title('Depth 30 km');
set(hTitle,'FontName','Helvetica', 'FontSize', 11, 'FontWeight','bold','units','inches');
P=get(hTitle,'Position');
P(2)=P(2)-.05;
set(hTitle, 'Position', P)
hColor=findobj(gcf,'tag','Colorbar');
C=get(hColor,'Position');
C{1}(1)=C{1}(1)-.01;
set(hColor(1),'Position',C{1});

subplot(3,2,6)
ray_plotmodel(X,Y,Z,w,'xy',40,'facet');
xlim([244 305])
ylim([212 292])
%xlim([240 310])
%ylim([200 300])
zlim([0 60])
shading interp
caxis([-4 4])
%caxis([-1.5 0.5])
set(gca,'Position',get(gca,'OuterPosition') - 1.3 * get(gca,'TightInset')...
        * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    
    % Uses prettier fonts than the defaults
set(gca,'FontName','Helvetica');
hTitle = title('Depth 40 km');
set(hTitle,'FontName','Helvetica', 'FontSize', 11, 'FontWeight','bold','units','inches');
P=get(hTitle,'Position');
P(2)=P(2)-.05;
set(hTitle, 'Position', P)
hColor=findobj(gcf,'tag','Colorbar');
C=get(hColor,'Position');
C{1}(1)=C{1}(1)-.01;
set(hColor(1),'Position',C{1});

set(gcf,'PaperPositionMode','auto')
print(gcf,'-painters','-depsc2','-r200','plot.eps')
%print(gcf,'-painters','-dpng','plot.png')