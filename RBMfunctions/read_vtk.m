function   [N_points,v,N_faces,f,c,c2,Tpsi,DTpsi,LabVent]=read_vtk(name2)
 
%Read vtk mesh with fibers (labels in points)
%N_points:number of vertex
%v: vertex coordinates
%N_faces: number of elements
%f: elements
%c: cell types
%

%name2='Fibers_0_0_20_kira32.vtk';

fileID = fopen(name2,'r');
%Read number of points
format_var='%s %f %s';
data = textscan(fileID,format_var,1,'HeaderLines',4,'Delimiter',' ');
N_points=data{1,2};
%save points
format_var=' %f %f %f';
data= textscan(fileID,format_var,N_points,'CollectOutput',1);
v=data{1,1};
%Read number of faces
format_var='%s %f %s';
data = textscan(fileID,format_var,1,'Delimiter',' ');
N_faces=data{1,2};
%save faces
format_var='%d %d %d %d %d';
data= textscan(fileID,format_var,N_faces,'CollectOutput',1);
f=data{1}(:,:);




%%%%%%%%%%%%%%
%cell types
format_var='%d';
dataempty=1;
 k=0;
  while dataempty==1 && k<10
    data = textscan(fileID,format_var,N_faces,'HeaderLines',k);
    dataempty=isempty(data{1,1});
   k=k+1;
  end
c=data{1,1};

%%
%read first type of data


 format_var='%s';
 dataempty=1;
 k=0;
  while dataempty==1 && k<10
     data = textscan(fileID,format_var,'HeaderLines',k);
    dataempty=isempty(data{1,1});
   k=k+1;
  end

 c2=data{1,1};
 %% header field
         format_header='%s';
         dataempty=1;
         k=0;
          while dataempty==1 && k<10
             data = textscan(fileID,format_header,1,'HeaderLines',k);
            dataempty=isempty(data{1,1});
           k=k+1;
          end
      
        format_header='%s';
         dataempty=1;
         k=1;
          while dataempty==1 && k<10
             data = textscan(fileID,format_header,1,'HeaderLines',k);
            dataempty=isempty(data{1,1});
           k=k+1;
          end

 %%
 format_var='%f %f %f';
 dataempty=1;
 k=1;
  while dataempty==1 && k<10
    data = textscan(fileID,format_var,N_points,'HeaderLines',k);
    dataempty=isempty(data{1,1});
   k=k+1;
  end
  DTpsi=cell2mat(data);
  %%
  format_var='%f';
  dataempty=1;
  k=0;
  while dataempty==1 && k<10
    data = textscan(fileID,format_var,N_points,'HeaderLines',k);
    dataempty=isempty(data{1,1});
   k=k+1;
  end
 Tpsi=data{1,1};
 
 dataempty=1;
  k=0;
  while dataempty==1 && k<10
    data = textscan(fileID,format_var,N_points,'HeaderLines',k);
    dataempty=isempty(data{1,1});
   k=k+1;  
  end  
 
 LabVent=data{1,1};
end
