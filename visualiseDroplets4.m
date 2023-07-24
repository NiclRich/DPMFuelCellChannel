%% Droplet visualisation model
% Take the x,y,z coordinates and make a 3D representation of the channel.
% patch the channel walls in 3D from Lx,Ly,Lz make Lz a maximum of 0.02
clear x y z

Lx = W;
Ly = H;
Lz = L;
gridRes = gridRes;

nx = gridRes+1;
ny = gridRes+1;
nz = gridRes.*(L/W)+1;

dx = Lx./nx;
dy = Ly./ny;
dz = Lz./nz;
% 
 x=0:dx:Lx;
 y=0:dy:Ly;
 z=0:dz:Lz;

%generate 3D grid for channel 
[X,Z,Y]=meshgrid(x,z,y);

%z=ones(size(X))
u=ones(size(X));
Wa=u.*0;
Regime=u.*0;


i=2:1:size(X,1)-1;
j=2:1:size(X,2)-1;
k=2:1:size(X,3)-1;

Circle=[D(:,3) D(:,4) D(:,5)];
Circle_Radius=D(:,7);
for kk=1:1:size(Circle,1)
    
distancemap=sqrt((X-Circle(kk,1)).^2 + (Y-Circle(kk,2)).^2 + (Z-Circle(kk,3)).^2);
mapfield=(Circle_Radius(kk,1)-distancemap);

for i=1:1:nz
    for j=1:1:nx
        for k=1:1:ny
             if mapfield(i,j,k)>0
             Wa(i,j,k)=1;
             Regime(i,j,k)=D(kk,6);
             end
        end
    end
end

end

%%

figure(1)


clf('reset')
s = Wa;
s = permute(s,[2 1 3]);

FV = isosurface(s); 
hold on
FV.vertices = FV.vertices.*[dx dy dz];

            p= patch('Faces',FV.faces,'Vertices',FV.vertices,'facecolor','b','linestyle','none');
            light('Position',[-1.5 1 1],'Style','local');
            axis('equal');
            view(3)
            set(p, 'DiffuseStrength',0.5, 'SpecularStrength',0.2, 'AmbientStrength',0.3);
            material("shiny");
            lighting flat

Lx = Lz;
Lz = Ly;

hold on
hold on
% make channel surfaces with corners
c1 = [0 0 0];
c2 = [0 0 Ly];
c3 = [0 Lz Ly];
c4 = [0 Lz 0];
wallPoints = [c1;c2;c3;c4];
view(3)
s3 = fill3(wallPoints(:,1)',wallPoints(:,2)',wallPoints(:,3)',wallPoints(:,3)','facecolor','white');
c1 = [0 0 0];
c2 = [Lx 0 0];
c3 = [Lx Lz 0];
c4 = [0 Lz 0];
wallPoints = [c1;c2;c3;c4];
s3 = fill3(wallPoints(:,1)',wallPoints(:,2)',wallPoints(:,3)',wallPoints(:,3)','facecolor','white');
c1 = [Lx 0 0];
c2 = [Lx 0 Ly];
c3 = [Lx Lz Ly];
c4 = [Lx Lz 0];
wallPoints = [c1;c2;c3;c4];
view(3)
s3 = fill3(wallPoints(:,1)',wallPoints(:,2)',wallPoints(:,3)',wallPoints(:,3)','facecolor','white','FaceAlpha',0.3);

c1 = [0 0 Ly];
c2 = [Lx 0 Ly];
c3 = [Lx Lz Ly];
c4 = [0 Lz Ly];
wallPoints = [c1;c2;c3;c4];
s3 = fill3(wallPoints(:,1)',wallPoints(:,2)',wallPoints(:,3)',wallPoints(:,3)','facecolor','white','FaceAlpha',0.3);

%daspect([1e-6 1e-6 1e-6]);
axis equal
%axis tight
%view(50,20)
camproj perspective
xlabel('x')
ylabel('y')
zlabel('z')
set(gcf, 'Position',  [500, 200, 1000, 1000]) %1400 %100
set(gca,'color','none')

hold on










