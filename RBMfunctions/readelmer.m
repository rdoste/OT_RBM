function   [N_points,v,N_faces,f,c,c2,Tpsi,DTpsi]=readelmer(name)

fileID = fopen(name,'r');
%Read number of points
format='%s %f %s';
data = textscan(fileID,format,1,'HeaderLines',4,'Delimiter',' ','CollectOutput',1);
N_points=data{1,2};
%save points
format=' %f %f %f';
data= textscan(fileID,format,N_points,'CollectOutput',1);
v=data{1,1};
%Read number of faces
format='%s %f %s';
data = textscan(fileID,format,1,'Delimiter',' ','CollectOutput',1);
N_faces=data{1,2};
%save faces
format='%d %d %d %d %d';
data= textscan(fileID,format,N_faces,'CollectOutput',1);
f=data{1}(:,:);
%%%%%%%%%%%%%%
format='%f';
data = textscan(fileID,format,N_faces,'HeaderLines',3);
c=data{1,1};
 format='%d';
 data = textscan(fileID,format,N_faces,'HeaderLines',5);
 c2=data{1,1};
  format='%f';
 data = textscan(fileID,format,N_points,'HeaderLines',4);
 Tpsi=data{1,1};
 format='%f';
 data = textscan(fileID,format,N_points*3,'HeaderLines',3);
 DTpsi1=data{1,1};
 format='%f';
 data = textscan(fileID,format,N_points*3,'HeaderLines',2);
 DTpsi2=data{1,1};
 format='%f';
 data = textscan(fileID,format,N_points*3,'HeaderLines',2);
 DTpsi3=data{1,1};
 
 DTpsi=[DTpsi1 DTpsi2 DTpsi3];
 
end