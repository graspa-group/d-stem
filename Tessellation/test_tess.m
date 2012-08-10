figure
subplot(1,2,1);
FV = sphere_tri('ico',3,1000);
lighting phong; shading interp; 
patch('vertices',FV.vertices,'faces',FV.faces,...
    'facecolor',[1 0 0],'edgecolor',[.2 .2 .6],'FaceAlpha',1);
axis off; camlight infinite; camproj('perspective');
axis equal
subplot(1,2,2);
FV2 = sphere_tri('ico',2,500);
lighting phong; shading interp; 
patch('vertices',FV2.vertices,'faces',FV2.faces,...
    'facecolor',[1 0 0],'edgecolor',[.2 .2 .6],'FaceAlpha',1);

axis off; camlight infinite; camproj('perspective');
axis equal