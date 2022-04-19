%   This script applies the OT_RBM to 3D meshes and obtains the fiber
%   orientation in each vertex of the mesh.
%
%     Copyright (C) 2019 Ruben Doste    Contact:ruben.doste@gmail.com
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.


clear
directory=pwd;
%Paths=====================================================================================
    if ismac
        disp('Platform not supported')
        elseif isunix  %add paths in Linux
        %elmer_path = '/usr/bin/';
        elmer_path='~/Programs/Elmer/install/bin';
                
             project_path = directory;
             directoryElmer=strcat(directory,'/Elmer/');
             project_path=strcat(project_path,'/Elmer/');
             directoryFunctions=strcat(directory,'/RBMfunctions/');
             addpath (genpath(directoryFunctions));    
        elseif ispc           %add paths in Windows
            elmer_path = '"C:\Program Files\Elmer 9.0-Release\bin"';
            directoryElmer=strcat(directory,'\Elmer\');
            project_path1 = directory;
            project_path=strcat(project_path1,'\Elmer\');
            directoryFunctions=strcat(directory,'\RBMfunctions\');
            addpath (genpath(directoryFunctions));

    end
   

%Variables=================================================================================================
      %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Remesh=0;
        name='Paloma28.vtk';
        namemesh=[name(1:end-3),'msh'];
      %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                          
        NewHeatCalculation=1;  %if 1, runs elmer calculations and neighbours calculation
      %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    
         ElementType=2;   %1--> Hexahedral
                          %2--> Tetrahedral
                          
                          
         largenumber=0;    %0--> no need interpolation for a large number of points
                           %1--> interpolation requiered after Elmer
                           
          name_vox='Paloma28.vtk'; %name of the hexahedral or the final tetrahedral mesh (only if largenumber is activated)
                           
             if ElementType==1
                 largenumber=1; %always interpolate in hexahedral
             end     
      %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    
        Trabeculations=1;    % Trabeculation detection
      %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    
       %Angle Definition
              Interpolation_in_septum=1; %1--> interpolate setpum points with continuity or disontinuity of fibers
                                         %0--> no interpolation
              Septum_angleRV=deg2rad(-0);%angle of the septal plane fibers respect to RV
              Discontinuity_angle=deg2rad(0);%angle between the RV septal fibers and LV septal fibers

                      AENDORV=90;
                      AEPIRV=-25;
                      AENDOLV=-60;
                      AEPILV=60;
                      AOTENDOLV=-90;
                      AOTENDORV=90;
                      AOTEPILV=0;
                      AOTEPIRV=0;
                  beta=1; % 1--> with transverse angle    0--> without angle    
                      beta_endoRV =deg2rad(180); 
                      beta_epiRV =deg2rad(180);  
                      beta_endoLV =deg2rad(180);
                      beta_epiLV=20;   %difference between epi/endo and the center of the LV
       %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                 
             ELVIRA=0;         %1--> save Elvira files 
       %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %% First Simulation   (only for calculating transmural gradient and determine RV and LV)
        %open vtk=================================================================================================

        [N_points_full,pto,N_faces_full,car,Fid]=read_labels2('labels.vtk');
               
        if Remesh==1
            %Remesh 
            Elementvol=1;
            [pt,p0,v0,t,idx]=surfinterior(pto,car+1);
            [node,elem,~]=surf2mesh(pto,car+1,[-100,0,0],[200,300,200],0.99,Elementvol,pt);
            opt.reratio=1.41;       
            elem=elem(:,1:end-1);        
            [node,elem,faces]=meshrefine(node,elem,[],opt);
            face=faces(:,1:end-1);
        else
            
             [N_points3d,node,N_faces3d,elem,~,~,~,~]=read_vtk(name);
             face=volface(elem(:,2:5));
             elem=elem(:,2:5)+1;
             face=face+1;
            
        end
       
%labelling

        if NewHeatCalculation==1                        
             centroid=meshcentroid(node,face);
             centroid_lab=meshcentroid(pto,car+1);
             TR=delaunayTriangulation(centroid_lab);
             [N,~]=nearestNeighbor(TR,centroid);
             label=Fid(N);
%                     %RV apex projection to epi                                                                                                                               -
%                     Epilab=find(Fid==2);                                                                                                                                      -
%                     TRepi=delaunayTriangulation(centroid_lab(Fid==2,:));                                                                                                      -
%                     [Napex,dapex]=nearestNeighbor(TRepi,centroid_lab(label==18,:));                                                                                            -
%                     label(Fid==18)=1;                                                                                                                                         -
%                     label(Epilab(Napex))=18;   
             centroid_tetra=meshcentroid(node,elem);
             [N_tetra,d_tetra]=nearestNeighbor(TR,centroid_tetra);
             label_tetra=Fid(N_tetra);
             label_tetra2=ones(size(label_tetra)); % we force only one physical body
             nnode=length(node);
             nelem=length(elem);
             carelem=size(elem,2);
             nface=length(face);
             save('labelsmesh','label','label_tetra');                 
        elseif NewHeatCalculation==0
             load('labelsmesh.mat');
            
        end
        
        
        %write msh file generation==============================================================================================
        if NewHeatCalculation==1
            mkdir Elmer;
            cd (directoryElmer)
            
            meshgenerator2_ot(namemesh,nnode,node,nelem,[zeros(nelem,1)+carelem elem],face,label,label_tetra2);

                %in case we don't have the pulmanary and aortic valve difined manually, we
                 %treat the rings as the valves
                 pulm=16;
                 aort=15;
                 if isempty(label==16)==1
                     pulm=9;
                 end
                 if isempty(label==15)==1
                     aort=10;
                 end

        %generate *.mesh files ==============================================================================================
                %delete old files if they exist
                if isequal(exist([project_path filesep 'mesh.header'],'file'),2)
                  try
                   delete([project_path filesep 'mesh.header']);
                    delete([project_path filesep 'mesh.nodes']);
                    delete([project_path filesep 'mesh.elements']);   
                    delete([project_path filesep 'mesh.boundary']);         
                  catch
                  end
                end

                mesh_command = [ elmer_path filesep 'ElmerGrid 14 2 ' namemesh ' -out ' project_path];  %run elmergrid.exe without arguments to get help on the command line parameters
                system(mesh_command);

        %generate *.sif file =================================================================================================
        %sif_generator(name,project_path,bodies,mat1,mat2,nboundaries,boundary1,boundary2,boundary3,temp1,temp2,temp3)
        %bodies can be 1 or 2   1-->99     2--> 56 and 78
        %boundaries goes from 1 to 3
        %mat1 and mat2 are the material heat conductivities (1 or 0)

                        %Transmural****************************************************************************************************************************
                        sif_generator_grad('Transmural',project_path,1,1,1,3,[2 12],[3 9 14 16 18],[1 10 13 15],0,1,-2);%000--> not defined
                        %Execute
                        simulation_command = ([ elmer_path filesep 'ElmerSolver '  project_path filesep 'case.sif']);   
                        system(simulation_command);
                        
                        %Longitudinal_bi***********************************************************************************************************************
                        sif_generator_grad('Longitudinal_bi',project_path,1,1,1,2,12,[9 14 16 10 13 15],000,1,0,000);%000--> not defined
                        %Execute
                        simulation_command = ([ elmer_path filesep 'ElmerSolver '  project_path filesep 'case.sif']);   
                        system(simulation_command);

                        %Transmural_bi ************************************************************************************************************************
                        sif_generator_grad('Transmural_bi',project_path,1,1,1,3,[2 12],[3 9 14 16 18],[1 10 13 15],0,1,1);%000--> not defined
                        %Execute
                        simulation_command = ([ elmer_path filesep 'ElmerSolver '  project_path filesep 'case.sif']);   
                        system(simulation_command);

           % Creation of a new surface: Septal wall                 

                      [N_points,v,N_faces1,f1,c,c2,Tphi,DTphi]=readelmer('transmural0001.vtk');  
                      f11=f1(:,1)==4;
                      f=f1(f11,:);
                      N_faces=length(f);
                            f2=f; %save original faces
                            v2=v;    %save original (future aha calculation)          
                            Tphi2=Tphi;
                            DTphi2=DTphi;
                      f1=f1(:,2:5)+1;

            [N2,dist]=nearestNeighbor(TR,v);         
            labelp=Fid(N2);

   %location of the septum surface creating two bodies (one for each
   %ventricle)
              body=Tphi(elem)./abs(Tphi(elem));
              body1=mode(body,2);
            
                  %in order to classify NaN, we look for the neighbors
                 TR_mesh=triangulation(double(f1),v2);
                  for iter=1:10
                      zeros_body=find(isnan(body1));
                      Neighbors_zeros = neighbors(TR_mesh,zeros_body);
                      Neighbors_val=body1(mode(Neighbors_zeros,2));
                      body1(zeros_body)=Neighbors_val;
                  end
              
              
              indxLV=find(body1<0);
              elemLV=elem(indxLV,:);
              surf_septL=volface(elemLV); 
              indxRV=find(body1>0);
              elemRV=elem(indxRV,:);
              surf_septR=volface(elemRV);
              points_sept=intersect(surf_septR,surf_septL);
              labelsept=labelp(points_sept);
              labelsept(labelsept==2 |labelsept==12)=22;
              points_sept2=points_sept(labelsept~=22);           
              points_RV=find(Tphi==1 | Tphi==-2 | Tphi ==0 ) ;
              points_RV_body=unique(surf_septR);
              points_Sept_RV=setdiff(points_RV_body,points_RV);         
              points_Sept=intersect(points_Sept_RV,points_sept2);
              label_sept=zeros(size(v));
              label_sept(points_Sept)=1;
              septal_surface_label=label_sept(surf_septR);
              septal_surface_index=find(prod(septal_surface_label,2));
              septal_surface=surf_septR(septal_surface_index,:);
              face2= [face;septal_surface];                      
              label2=[label;repmat(20,length(septal_surface),1)]; %          
              label_tetra_msh=ones(size(label_tetra)).*10;
              label_tetra_msh(body1<0)=11;                  

          meshgenerator2_ot(namemesh,nnode,node,nelem,[zeros(nelem,1)+carelem elem],face2,label2,label_tetra_msh); 

                 %generate *.mesh files for longitudinal gradients=========================================================================================
                            %delete old files if they exist
                            if isequal(exist([project_path filesep 'mesh.header'],'file'),2)
                              try
                               delete([project_path filesep 'mesh.header']);
                                delete([project_path filesep 'mesh.nodes']);
                                delete([project_path filesep 'mesh.elements']);   
                                delete([project_path filesep 'mesh.boundary']);         
                              catch
                              end
                            end

                            mesh_command = [ elmer_path filesep 'ElmerGrid 14 2 ' namemesh ' -out ' project_path];  %run elmergrid.exe without arguments to get help on the command line parameters
                            system(mesh_command);

        %Long_RV*******************************************************************************************************************************
                            sif_generator_flux('Longrv',project_path,2,0,1,2,18,[9 14 16],000,1,0,000);%000--> not defined
                            %Execute
                            simulation_command = [ elmer_path filesep 'ElmerSolver '  project_path filesep 'case.sif'];   
                            system(simulation_command);

                            %Long_LV*******************************************************************************************************************************
                            sif_generator_flux('Longlv',project_path,2,1,0,2,12,[10 13 15],000,1,0,000);%000--> not defined
                            %Execute
                            simulation_command = [ elmer_path filesep 'ElmerSolver '  project_path filesep 'case.sif'];   
                            system(simulation_command);

                            %2valveRV*******************************************************************************************************************************
                            sif_generator_flux('2valverv',project_path,2,0,1,3,18,14,9,1,1,0);%000--> not defined
                            %Execute
                            simulation_command = [ elmer_path filesep 'ElmerSolver '  project_path filesep 'case.sif'];   
                            system(simulation_command);

                            %2valveRV*******************************************************************************************************************************
                            sif_generator_flux('2valvelv',project_path,2,1,0,3,12,13,10,1,1,0);%000--> not defined
                            %Execute
                            simulation_command = [ elmer_path filesep 'ElmerSolver '  project_path filesep 'case.sif'];   
                            system(simulation_command);

                            %Septal*******************************************************************************************************************************
        %                     sif_generator_flux('Septum',project_path,2,1,1,3,[2 12],[3 9 14 pulm 18 1 10 13 aort],20,1,1,0);%000--> not defined
        %                     %Execute
        %                     simulation_command = [ elmer_path filesep 'ElmerSolver '  project_path filesep 'case.sif'];   
        %                     system(simulation_command);

                            sif_generator_flux('Septum',project_path,2,1,1,2,[3 9 14 16 18 1 10 13 15],20,000,0,1,000);%000--> not defined
                            %Execute
                            simulation_command = [ elmer_path filesep 'ElmerSolver '  project_path filesep 'case.sif'];   
                            system(simulation_command);

         %% Fiber Generation 

        elseif  NewHeatCalculation==0
            cd Elmer    
               [N_points,v,N_faces1,f1,c,c2,Tphi,DTphi]=readelmer('transmural0001.vtk');  

                      f11=f1(:,1)==4;
                      f=f1(f11,:);
                      N_faces=length(f);
                            f2=f; %save original faces
                            v2=v;    %save original (future aha calculation)          
                            Tphi2=Tphi;
                            DTphi2=DTphi;
                      f1=f1(:,2:5)+1;                
        end          

         %[N_points,v,N_faces1,f1,c,c2,Tpsi,DTpsi]=readelmer('longitudinal_base0001.vtk');
         %[N_points,v,N_faces1,f1,c,c2,Tphi,DTphi]=readelmer('transmural0001.vtk');
         [~,~,~,~,~,~,Tval,~]=readelmer('2valverv0001.vtk');
         [~,~,~,~,~,~,Tval_l,~]=readelmer('2valvelv0001.vtk');%interpolare con este valor
         [~,~,~,~,~,~,T_sept,DT_sept]=readelmer('septum0001.vtk');%interpolacion septo (fiber discontinuity)
         [~,~,~,~,~,~,TpsiR,DTpsiR]=readelmer('longrv0001.vtk');%Longitudinal gradient (we find the closest point)
         [~,~,~,~,~,~,TpsiL,DTpsiL]=readelmer('longlv0001.vtk');%Longitudinal gradient (we find the closest point)
         [~,~,~,~,~,~,Tphi_bi,DTphi_bi]=readelmer('transmural_bi0001.vtk');  
         [~,~,~,~,~,~,Tpsi_bi,DTpsi_bi]=readelmer('longitudinal_bi0001.vtk');  
       
    %normalize gradients and change sign in flux

    DTphi=double(normr(DTphi));
    DTpsiL1=DTpsiL;
    DTpsiR1=DTpsiR;%old ones, they will be useful later
    DTpsiL=double(normr(DTpsiL).*-1);
    DTpsiR=double(normr(DTpsiR).*-1);
    Tval(Tval<0)=0;
    Tval_l(Tval_l<0)=0; 
    Tval2=Tval; %future aha calculation
     
   %Trabeculations
          if Trabeculations==1
              Magrad= sqrt(sum(DTphi2.^2,2));    
          end
          
  %%%%%% We define LV or RV tags (epi-endo)%%%%%%%%%%%%%%%%%%%%%%%%%
      Ventricle=Tphi; 
      Ventricle(Ventricle<0 & Ventricle>-2)=-1; % mio LV =-1
      Ventricle(Ventricle>0 & Ventricle<1)=2;  % mio RV =2 
      gradientdot=dot(DTphi,DTphi_bi,2);%the different direction of the gradients will determine LV vs RV
      Ventricle(Tphi==0 & gradientdot>0)=3;
      Ventricle(Tphi==0 & gradientdot<0)=-3;
 
 cd ..
   %% Interpolation for large numbers
                 if largenumber==1
                      if ElementType==1                 
                        [N_points_vox,v_vox,N_faces_vox,f_vox]=readvoxel(name_vox);                       
                      else
                             [N_points_vox,v_vox,N_faces_vox,f_vox,~,~,~,~]=read_vtk(name);
                            
                      end
                          
                        node_vox=v_vox;
                    
                        elem_vox=f_vox;
                        
                        %face_iris=volface(elem_iris(:,2:5));
                        %f=elem_iris;
                        %elem_iris=elem_iris(:,2:5)+1;
                        %face_iris=face_iris+1;

                        % Laplace equation interpolation
                        
                        T=triangulation(double(f2(:,2:5)+1),double(v2));
                        tic
                        [tri,bar]=pointLocation(T,node_vox(:,:));
                        toc
                        %substitute NaN by nearestNeighbor
                        [row, col] = find(isnan(tri));
                        pointN=(nearestNeighbor(T,node_vox(row,:)));
                        [n1,n2]=ismember(pointN,f2(:,2:end)+1);
                        [row3, ~] = ind2sub(size(f2(:,2:end)+1), n2);
                        tri(row)=row3;
                                               
                       bar(row,:)=cartesianToBarycentric(T,row3,v2(pointN,:));
                    
                        points_Tetra=f2(tri,2:end)+1;
                        
                        % Values interpolation
                        TpsiRVpointsTetra=TpsiR(points_Tetra);
                        TpsiR=dot(bar,TpsiRVpointsTetra,2);
                        
                        TpsiLVpointsTetra=TpsiL(points_Tetra);
                        TpsiL=dot(bar,TpsiLVpointsTetra,2);

                        TphipointsTetra=Tphi(points_Tetra);
                        Tphi=dot(bar,TphipointsTetra,2);
                        
                        TphibipointsTetra=Tphi_bi(points_Tetra);
                        Tphi_bi=dot(bar,TphibipointsTetra,2);

                        TpsibipointsTetra=Tpsi_bi(points_Tetra);
                        Tpsi_bi=dot(bar,TpsibipointsTetra,2);


                        TpsiLpointsTetra=T_sept(points_Tetra);
                        T_sept=dot(bar,TpsiLpointsTetra,2);

                        TvalpointsTetra=Tval(points_Tetra);
                        Tval=dot(bar,TvalpointsTetra,2);
                        
                        TvallpointsTetra=Tval_l(points_Tetra);
                        Tval_l=dot(bar,TvallpointsTetra,2);

                        %  Interpolacion de gradientes
                         %DTpsiR
                        DTpsiR=gradientinterpol_noGPU(DTpsiR1,length(node_vox),points_Tetra,bar);
                        DTpsiL=gradientinterpol_noGPU(DTpsiL1,length(node_vox),points_Tetra,bar);
                        DTphi=gradientinterpol_noGPU(DTphi,length(node_vox),points_Tetra,bar);
                        DTphi_bi=gradientinterpol_noGPU(DTphi_bi,length(node_vox),points_Tetra,bar);
                        DTpsi_bi=gradientinterpol_noGPU(DTpsi_bi,length(node_vox),points_Tetra,bar);
                        %change points of calculations
                        v=node_vox;
                        N_points=length(node_vox);
                        f=elem_vox;
                        f1=f(2:end)-1;
                        N_faces=length(elem_vox);
                        name=name_vox;
                         
                
                %Ventricle Tag interpolation
                        
                          Ventricle2=int8(Ventricle(points_Tetra));
                          Ventricle_raw=round(dot(bar,double(Ventricle2),2));%know lV vs RV 
                          %rounded value
                          Ventricle3=Ventricle2;%auxiliar
                          %positive values
                          Ventricle3(Ventricle2==-2)=1;
                          Ventricle3(Ventricle2==-1)=2;
                          Ventricle3(Ventricle2==-3)=3;
                          %Ventricle values calculation
                          Ventricle=round(dot(bar,double(Ventricle3),2));
                          Ventricle(Ventricle_raw<0)=Ventricle(Ventricle_raw<0).*-1;            
                          Ventricle3=int8(Ventricle);%auxiliar
                          Ventricle(Ventricle3==-1)=-2;
                          Ventricle(Ventricle3==-2)=-1;
                   clear -regexp pointsTetra$;                        
                        
                 end
                
  %% Axis
    % e1 longitudinal+vertical(apex 2 mitral and triculpid valves)
        DTpsi_result=normr(DTpsiL);  %2ValvesDirection+PVDirection
        indR=find(Ventricle==1 | Ventricle==2 | Ventricle==3); %RV
        DTpsi_result(indR,:)=normr(DTpsiR(indR,:));
        e1=normr(DTpsi_result);
        e1=normr(e1);

    %%%%%% DTphi_final is an interpolotaion between the two transmural directions
    %intentar bislerp (to do)
        DTphi3=DTphi;
        DTphi(Ventricle<0,:)=DTphi(Ventricle<0,:).*-1;
        DTphi_final=normr(Tphi_bi.*DTphi+(1-Tphi_bi).*DTphi_bi);
        DTphi_final(Ventricle<0,:)=DTphi_final(Ventricle<0,:).*-1;
        proj=DTphi_final-bsxfun(@times,dot(e1',DTphi_final')',e1);
    %%%%%%%%%%%%

    %proj=DTphi-bsxfun(@times,dot(e1',DTphi')',e1);
        e2=normr(proj);
        e2=normr(e2);
        e0=cross(e1',e2')';
        e0=normr(e0);
        e0=normr(e0); %uso normr 2 veces para normalizar los vectores que son (0,0,0)-->(1,1,1)-->Normalizados (not important)
        e1_bi=normr(DTpsi_bi);
        proj_bi=DTphi_final-bsxfun(@times,dot(e1_bi',DTphi_final')',e1_bi);
        e2_bi=normr(proj_bi);
        e0_bi=cross(e1_bi',e2_bi')';
        e0_bi=normr(e0_bi);

    %Defino transmuralidad (0-epi   1-endo)
    d=abs(abs(Tphi)-1); %invierto el orden en d. Asi d=0 endo d=1 epi
    d2=abs(abs(Tphi)/2-1); %for LV
    d_bi=Tphi_bi;

%%  alpha angle 

    alpha_wall=zeros(size(Ventricle));
    alpha_wall_epi=zeros(size(Ventricle));
    alpha_wall_endo=zeros(size(Ventricle));     
    
    %Orientacion Epicardio               
          alpha_wall(Ventricle==-3)=deg2rad(AEPILV).*Tval_l(Ventricle == -3) +deg2rad(AOTEPILV) .*(1-Tval_l(Ventricle==-3));   %OT interpolation
          alpha_wall(Ventricle==3 )=deg2rad(AEPIRV).*Tval (Ventricle == 3)   +deg2rad(AOTEPIRV) .*(1-Tval(Ventricle==3));   %OT interpolation
          
    %Orientacion Endocardio 
          alpha_wall(Ventricle==-2)=deg2rad(AENDOLV).*Tval_l(Ventricle==-2)  +deg2rad(AOTENDOLV) .*(1-Tval_l(Ventricle==-2));  %OT interpolation
          alpha_wall(Ventricle==1 )=deg2rad(AENDORV).*Tval(Ventricle==1)     +deg2rad(AOTENDORV) .*(1-Tval(Ventricle==1));     %OT interpolation 
        
    %Orientacion Miocardio
   
        %RV
   
       indx=find(Ventricle==2);  %interpolacion RV
            alpha_wall_epi(indx)=deg2rad(AEPIRV).*Tval (indx)   +deg2rad(AOTEPIRV) .*(1-Tval(indx));  %OT interpolation
            alpha_wall_epi(indx)=alpha_wall_epi(indx).*(1-T_sept(indx))+T_sept(indx).*(Septum_angleRV+alpha_wall_epi(indx)); %Septum interpolation
            alpha_wall_endo(indx)=deg2rad(AENDORV).*Tval(indx)  +deg2rad(AOTENDORV).*(1-Tval(indx));  %OT interpolation 
       
        alpha_wall(indx)=alpha_wall_endo(indx).*(1-d(indx))+alpha_wall_epi(indx).*(d(indx));   %linear transmural interpolation       
       
        %LV
   
       indx=find(Ventricle==-1); % interpolacion LV
            alpha_wall_epi(indx)=deg2rad(AEPILV).*Tval_l (indx)   +deg2rad(AOTEPILV) .*(1-Tval_l(indx));  %OT interpolation
         %interpolation of septum comes later         
            alpha_wall_endo(indx)=deg2rad(AENDOLV).*Tval_l(indx)  +deg2rad(AOTENDOLV).*(1-Tval_l(indx));  %OT interpolation
        
        alpha_wall(indx)=alpha_wall_endo(indx).*(1-d2(indx))+alpha_wall_epi(indx).*(d2(indx));   %linear transmural interpolation
          
%% 
        if Interpolation_in_septum==1
    %    Correction of the septum points
    %     this step is neccesary in order to achieve a continuity in the septum.
    %     We can chooose to modify the rv septal points or lv septal points.
    %     Here we modify the LV points (more points)
    %     For calculating this "angle of correction", we have to calculate
    %     individually for each septal point, and compare the e1 vector of the
    %     LV with the e1 vector of the RV 
            points_sept_all=find(T_sept>0.15);
            points_sept_RV=intersect(find(Ventricle==2),points_sept_all);
            points_sept_LV=intersect(find(Ventricle==-1 | Ventricle==-3),points_sept_all);
            TR_sept=delaunayTriangulation(v(points_sept_RV,:));
            [N_sept,~]=nearestNeighbor(TR_sept,v(points_sept_LV,:));
            e1_rvsept=e1(points_sept_RV,:);
            e1_rvsept_ref=e1_rvsept(N_sept,:);
            e1_lvsept=e1(points_sept_LV,:);
            prod_vect=cross(e1_rvsept_ref,e1_lvsept,2);
            angle_correct=atan2(sqrt(sum(prod_vect.^2,2)),dot( e1_rvsept_ref,e1_lvsept,2));%angle between both e1 vectors

            Septum_angleLV=angle_correct+Discontinuity_angle;        
            
            indx=points_sept_LV; %  LV interpolation
                alpha_wall_epi(indx)=deg2rad(AEPILV).*Tval_l (indx)   +deg2rad(AOTEPILV) .*(1-Tval_l(indx));  %OT interpolation
                alpha_wall_epi(indx)=alpha_wall_epi(indx).*(1-T_sept(indx))+T_sept(indx).*(Septum_angleLV+Septum_angleRV); %Septum interpolation
                alpha_wall_endo(indx)=deg2rad(AENDOLV).*Tval_l(indx)  +deg2rad(AOTENDOLV).*(1-Tval_l(indx));  %OT interpolation
                alpha_wall(indx)=alpha_wall_endo(indx).*(1-d2(indx))+alpha_wall_epi(indx).*(d2(indx)); 
        end
          
                        
             %% beta angle
                 
          beta_epiLVext =deg2rad(180-beta_epiLV);
          beta_epiLVint =deg2rad(180+beta_epiLV);
          beta_wall=beta_endoRV*(1-d)+beta_epiRV*(d);   %all heart (RV included)
%         beta_wall(Ventricle<0)=-beta_endo*(1-d2(Ventricle<0))-beta_epi*d2(Ventricle<0); 
          %d22=abs(abs(d)-1);    %endo0   mid 1    epi 0
          %only for LV
          beta_wall(Tphi<=-1)=beta_endoLV*(1-d(Tphi<=-1))+ beta_epiLVint*(d(Tphi<=-1));
          beta_wall(Tphi>-1 & Tphi<0 | Ventricle==-3 )=beta_endoLV*(1-d(Tphi>-1 & Tphi<0 | Ventricle==-3 ))+beta_epiLVext*(d(Tphi>-1 & Tphi<0 | Ventricle==-3 ));
          
          beta_wall=beta_wall.*abs(T_sept-1)+deg2rad(180).*T_sept; %interpolacion septo
          beta_wall(Ventricle>0)=beta_wall(Ventricle>0).*Tval(Ventricle>0)+deg2rad(180).*(1-Tval(Ventricle>0));%interpolacion OT
          beta_wall(Ventricle<0)=beta_wall(Ventricle<0).*Tval_l(Ventricle<0)+deg2rad(180).*(1-Tval_l(Ventricle<0));%interpolacion OT
          
          beta_wall=beta_wall*beta;  

              %Trabeculations
          if Trabeculations==1
              if largenumber==1
                  Magrad=Magrad(points_Tetra);
                  Magrad2=dot(bar, Magrad,2);
              else
                  Magrad2=Magrad;
              end
              
              Trabs_raw1=find(Magrad2<0.1);
              Trabs_raw2=find(Tphi>0.8 | Tphi<-1.8);
              Trabs=intersect(Trabs_raw1,Trabs_raw2); 
              Trabs_raw3=find(TpsiR>0.05 | TpsiL>0.05); %not in the rings             
              Trabs=intersect(Trabs,Trabs_raw3);% trabecular points
%               Trabs_label=zeros(size(Tphi));
%               Trabs_label(Trabs)=1;             
              MagTrab=abs(Magrad2(Trabs))/max(Magrad2(Trabs));
              alpha_wall(Trabs)=deg2rad(90)*(1-MagTrab)+alpha_wall(Trabs).*MagTrab;
              beta_wall(Trabs)=0;
          end
          
          %% Rotation              
  %creation of 3x3xN  matrix with vectors e0 e1 e2 
      Q=cat(3,e0,e1,e2); 
      Q=permute(Q,[3 2 1]);
     
    %ration matrix creation
   %rotacion over e2
      axis=e2;
      theta=alpha_wall;
      Rot_axis=[cos(theta)+(axis(:,1).^2).*(1-cos(theta)) axis(:,1).*axis(:,2).*(1-cos(theta))-axis(:,3).*sin(theta) axis(:,1).*axis(:,3).*(1-cos(theta))+axis(:,2).*sin(theta);
             axis(:,2).*axis(:,1).*(1-cos(theta))+axis(:,3).*sin(theta) cos(theta)+(axis(:,2).^2).*(1-cos(theta)) axis(:,2).*axis(:,3).*(1-cos(theta))-axis(:,1).*sin(theta);
             axis(:,3).*axis(:,1).*(1-cos(theta))-axis(:,2).*sin(theta) axis(:,3).*axis(:,2).*(1-cos(theta))+axis(:,1).*sin(theta) cos(theta)+(axis(:,3).^2).*(1-cos(theta))];

        R=reshape(Rot_axis',3,N_points,3);
     R=permute(R,[3 1 2]);
     
      QX=zeros(size(R));
      for indx=1:length(alpha_wall)
          QX(:,:,indx)= mtimes(Q(:,:,indx), R(:,:,indx));
      end
  % second rotation over ef
      axis2=QX(2,:,:);
      theta2=beta_wall;
      axis=reshape(axis2,3,N_points)';
      Rot_axis2=[cos(theta2)+(axis(:,1).^2).*(1-cos(theta2)) axis(:,1).*axis(:,2).*(1-cos(theta2))-axis(:,3).*sin(theta2) axis(:,1).*axis(:,3).*(1-cos(theta2))+axis(:,2).*sin(theta2);
             axis(:,2).*axis(:,1).*(1-cos(theta2))+axis(:,3).*sin(theta2) cos(theta2)+(axis(:,2).^2).*(1-cos(theta2)) axis(:,2).*axis(:,3).*(1-cos(theta2))-axis(:,1).*sin(theta2);
             axis(:,3).*axis(:,1).*(1-cos(theta2))-axis(:,2).*sin(theta2) axis(:,3).*axis(:,2).*(1-cos(theta2))+axis(:,1).*sin(theta2) cos(theta2)+(axis(:,3).^2).*(1-cos(theta2))];

        
    R2=reshape(Rot_axis2',3,N_points,3);
    R2=permute(R2,[3 1 2]);   
  
    QX2=zeros(size(R));

    for indx=1:length(alpha_wall)
        QX2(:,:,indx)= mtimes(QX(:,:,indx), R2(:,:,indx));
    end

    FX=QX2(1,:,:);
    F=reshape(FX,3,N_points)';

clear -regexp ^Q;
clear -regexp ^R;          
        %% Final "make-up" (valves, apex presentation...)  
 
    % Select if you want to show the angle or the ot
    alpha_wall2=rad2deg(alpha_wall);
    beta_wall2=rad2deg(beta_wall);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %OT definition (for conductivity purposes)  
                    OutflowT=ones(size(Ventricle));
                    OutflowT((Tval_l<0.05 & Ventricle < 0) | (Tval_l<0.05 & Ventricle==-3)  )=2; %AorticOT
                    OutflowT((Tval<0.05 & Ventricle >0) | (Tval<0.05 & Ventricle==3) )=3;       %Pulmonary OT
                    OutflowT((Tval_l==0 & Ventricle  < 0) | (Tval_l==0 & Ventricle==-3)  )=4; %Aortic Valve             
                    OutflowT((Tval==0 & Ventricle  > 0) | (Tval==0 & Ventricle==3)  )=5; %Pulmanoary Valve
                    OutflowT(((Tval_l==1 & TpsiL~=1) & Ventricle  < 0) | ((Tval_l==1 & TpsiL~=1) & Ventricle==-3)  )=6; %Mitral
                    OutflowT(((Tval==1 & TpsiR~=1) & Ventricle  > 0) | ((Tval==1 & TpsiR~=1) & Ventricle==3)  )=7;      %Tricuspid

                if largenumber==1
                    OutflowT=ones(size(Ventricle));
                    OutflowT((Tval_l<0.05 & Ventricle < 0) | (Tval_l<0.05 & Ventricle==-3)  )=2; %AorticOT
                    OutflowT((Tval<0.05 & Ventricle >0) | (Tval<0.05 & Ventricle==3) )=3;       %Pulmonary OT
                    OutflowT((abs(Tval_l)<0.001 & Ventricle  < 0) | (Tval_l==0 & Ventricle==-3)  )=4; %Aortic Valve             
                    OutflowT((abs(Tval)<0.001 & Ventricle  > 0) | (Tval==0 & Ventricle==3)  )=5; %Pulmanoary Valve
                    OutflowT(((abs(Tval_l)>0.997  & TpsiL<0.997) & Ventricle  < 0) | ((abs(Tval_l)>0.997  &  TpsiL<0.997) & Ventricle==-3)  )=6;
                    OutflowT(((abs(Tval)>0.997  & TpsiR<0.997) & Ventricle  > 0) | ((abs(Tval)>0.997  &  TpsiR<0.997) & Ventricle==3)  )=7;

                    %smooth OutflowT (only tetra)
%                      Outflow_elem=OutflowT(f1);
%                      Outflow_elem=max(Outflow_elem,[],2);
%                      for Outflowlabel=2:7
%                      Outflow_elem=smoothlabel(Outflow_elem,Outflowlabel,1,'faces',f1);
%                      end
% 
%                      [~,uuu]=unique(f1);
%                      [row4, col4] = ind2sub(size(f1), uuu);
%                      OutflowT=Outflow_elem(row4);
%                      OutflowT(OutflowT>1)=2;
                end  

     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 


% transmurality

d3=abs(abs(Tphi)-1); 
d3(Tphi<0)=abs(abs(Tphi(Tphi<0))/2-1);

%% 33% epi, endo mio
   [~,indpc1,~]=unique(f(:,2:end));  %labels tetra --> labels points
   if size(f1,2)==4
    
    label_tetra3=repmat(label_tetra,1,4);
    labelpc=label_tetra3(indpc1);

    Ventricle3=ones(size(d3))*2;
    Ventricle3(d3<0.17)=1;
    Ventricle3(d3>0.60 & (labelpc==2 | labelpc==12))=3;
   elseif size(f1,2)==8
       
       Ventricle3=ones(size(Tphi_bi))*2;
       Ventricle3(Tphi_bi<0.17)=1;
       Ventricle3(Tphi_bi>0.83)=3;
    
   end


 
%%%%%

%diferent labels assignation to final mesh 
% possible labels: d3,Ventricle3, OutflowT, alpha_wall2, beta_wall2...

    label1=Ventricle;
    label3=Ventricle3;

    if Trabeculations==1    
        d4=ones(size(d3));
        d4(Trabs)=0;
        label2=d4;
    else
       label2=d3;
    end



 %save('Costfunction','Ventricle','e0','e1','e2','Tphi','Tval','Tval_l','T_sept','d2','v','d','aha','N_points');   

    string = num2str(round(rad2deg(Discontinuity_angle)));
    string = strcat(string,'_');
    string = strcat(string,num2str(round(rad2deg(Septum_angleRV))));
    string = strcat(string,'_');
    string = strcat(string,num2str(beta_epiLV));
    string = strcat(string,'_');
        if Interpolation_in_septum==0
        string='disc_';
        end
    name3=strcat('Fibers_',string); 
    name3=strcat(name3,name); 

 %write vtk file
 if ElementType==1
    write_vtk_rbm_hex(name3,N_points,v,N_faces,f,label1,label2,label3,F);    %write mesh with fibers
 elseif ElementType==2
     write_vtk_rbm(name3,N_points,v,N_faces,f,label1,label2,label3,F);
 end
                 






