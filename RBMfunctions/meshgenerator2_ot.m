%This function creates a tethraedric mesh (.msh file) wich serves as input file for
%ELMER FEM
function   meshgenerator2_ot(name,N_points,points,N_faces,faces,surface,label,label_tetra)

%name       output name
%boundary   differents element for elmer boundary condition (integer)
%class      type of element (tetragonal--4 ,hexahedral--5, triangles2D--2)



class2=faces(:,1);
class3=class2;
class3(find(class2==8))=5;
class3(find(class2==4))=4;
class3(find(class2==3))=2;

faces=faces(:,2:end);
s=size(faces,2);
ss=repmat(' %d', 1, s);

body=class3*0+99;
 body(find(label_tetra==5 | label_tetra==6 ))=56;
 body(find(label_tetra==7 | label_tetra==8 ))=78;
 
body(find(label_tetra==10))=78;
body(find(label_tetra==11))=56;
%% añadimos las caras
Nsurface=size(surface,1);
Classurface=size(surface,2);
Classurface=zeros(size(surface(:,1)))+Classurface;
Classurface(find(Classurface==8))=5;
Classurface(find(Classurface==4))=4;
Classurface(find(Classurface==3))=2;

s=size(surface,2);
ss2=repmat(' %d', 1, s);
%%
fileID = fopen(name, 'w');
fprintf(fileID,'$MeshFormat\n' );
fprintf(fileID,' %d %d %d \n',[2 0 8]);
fprintf(fileID,'$EndMeshFormat \n' );
fprintf(fileID,'$Nodes\n');
fprintf(fileID,'%d \n',N_points);
pointsm=[(1:N_points)',points];
fprintf(fileID,'%d %8.3f %8.3f %8.3f \n',pointsm');
fprintf(fileID,'\n');
fprintf(fileID,'$EndNodes\n');
fprintf(fileID,'$Elements\n');
fprintf(fileID,'%d \n',N_faces+Nsurface);
facesm=[(1:N_faces)',class3,repmat(2,(N_faces),1),body,class3,faces];
fprintf(fileID,strcat('%d %d %d %d %d ', ss ,' \n'),facesm');
surfacesm=[(N_faces+1:N_faces+Nsurface)',Classurface,repmat(2,(Nsurface),1),label,label,surface];
fprintf(fileID,strcat('%d %d %d %d %d ', ss2 ,' \n'),surfacesm');
fprintf(fileID,'$EndElements\n');
fclose(fileID);
end