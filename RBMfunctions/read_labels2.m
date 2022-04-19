function   [N_points_full,pto,N_faces_full,car,Fid]=read_labels2(name)

fileID = fopen(name,'r');
%Read number of points
format='%s %f %s';
data = textscan(fileID,format,1,'HeaderLines',4,'Delimiter',' ');
N_points_full=data{1,2};
%save points
format=' %f %f %f';
data= textscan(fileID,format,N_points_full,'CollectOutput',1);
pto=data{1,1};
%Read number of faces
format='%s %f %f';
data = textscan(fileID,format,1,'Delimiter',' ');

N_faces_full=data{1,2}(1);
if isnan(N_faces_full)==1
   data = textscan(fileID,format,1,'Delimiter',' ');
end
    
connectivity=data{1,3}/N_faces_full;
%save faces
format=repmat('%d ',1,connectivity);
data= textscan(fileID,format,N_faces_full,'CollectOutput',1);
car=data{1}(:,2:end);


format='%f';
data= textscan(fileID,format,N_faces_full,'HeaderLines',4);
k=1;
while isempty(data{1,1})==1 && k<50
   data= textscan(fileID,format,N_faces_full,'HeaderLines',k);
   k=k+1;
end
    
Fid=data{1,1};

fclose(fileID);

end