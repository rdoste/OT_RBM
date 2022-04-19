function   [N_points_vox,v_vox,N_faces_vox,f_vox]=readvoxel(name)
   

    fileID = fopen(name,'r');
    format='%s %f %s';
    data = textscan(fileID,format,1,'HeaderLines',4,'Delimiter',' ','CollectOutput',1);
    N_points_vox=data{1,2};
    %save points
    format=' %f %f %f';
    data= textscan(fileID,format,N_points_vox,'CollectOutput',1);
    v_vox=data{1,1};
    %Read number of faces
    format='%s %f %s';
    data = textscan(fileID,format,1,'Delimiter',' ','CollectOutput',1);
    N_faces_vox=data{1,2};
    %save faces
    format='%d %d %d %d %d %d %d %d %d';
    data= textscan(fileID,format,N_faces_vox,'CollectOutput',1);
    f_vox=data{1}(:,2:9);
    clear data
end