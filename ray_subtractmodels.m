function r = ray_subtractmodels(starting,final)

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

for i=1:length(r.modelijk(:,4))
    r.modelijk(i,4)=(starting.modelijk(i,4)-final.modelijk(i,4))*100/starting.modelijk(i,4);
    r.modelxyz(i,4)=(starting.modelxyz(i,4)-final.modelxyz(i,4))*100/starting.modelxyz(i,4);
    r.modelijk(i,6)=(starting.modelijk(i,6)-final.modelijk(i,6));
    r.modelxyz(i,6)=(starting.modelxyz(i,6)-final.modelxyz(i,6));
end