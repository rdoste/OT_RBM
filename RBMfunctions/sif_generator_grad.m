  function  sif_generator_grad(name,project_path,bodies,mat1,mat2,nboundaries,boundary1,boundary2,boundary3,temp1,temp2,temp3)
  
  %bodies can be 1 or 2   1-->99     2--> 56 and 78
  %boundaries goes from 1 to 3
  %mat1 and mat2 are the material heat conductivities (1 or 0)
  
  lb = sprintf('\n');
  
        %generate *.sif file =================================================================================================
        %define path
        sif_path =  [project_path filesep 'case.sif'];

        %delete old file if it exists
        if isequal(exist(sif_path,'file'),2)
          try
           delete(sif_path);
          catch
          end
        end

           %open file
        fileid = fopen(sif_path,'w');

        %write file
        fwrite(fileid,['Header' lb]);
        fwrite(fileid,[' CHECK KEYWORDS Warn' lb]);
        %~ fwrite(fileid,[' Mesh DB "'  project_path  '"' lb]);
        %~ fwrite(fileid,[' Include Path "' project_path '"' lb]);
        %~ fwrite(fileid,[' Results Directory "' project_path '"' lb]);
        fwrite(fileid,[' Mesh DB "." "."' lb]);
        fwrite(fileid,[' Include Path ""' lb]);
        fwrite(fileid,[' Results Directory ""' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);

        fwrite(fileid,['Simulation' lb]);
        fwrite(fileid,[' Max Output Level = 5' lb]);
        fwrite(fileid,[' Coordinate System = Cartesian' lb]);
        fwrite(fileid,[' Coordinate Mapping(3) = 1 2 3' lb]);
        fwrite(fileid,[' Simulation Type = Steady state' lb]);
        fwrite(fileid,[' Steady State Max Iterations = 1' lb]);
        fwrite(fileid,[' Output Intervals = 1' lb]);
        fwrite(fileid,[' Timestepping Method = BDF' lb]);
        fwrite(fileid,[' BDF Order = 1' lb]);
        fwrite(fileid,[' Solver Input File = case.sif' lb]);
        fwrite(fileid,[' Post File = case.ep' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);

        fwrite(fileid,['Constants' lb]);
        fwrite(fileid,[' Gravity(4) = 0 -1 0 9.82' lb]);
        fwrite(fileid,[' Stefan Boltzmann = 5.67e-08' lb]);
        fwrite(fileid,[' Permittivity of Vacuum = 8.8542e-12' lb]);
        fwrite(fileid,[' Boltzmann Constant = 1.3807e-23' lb]);
        fwrite(fileid,[' Unit Charge = 1.602e-19' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);
        
        
        if bodies==1

        fwrite(fileid,['Body 1' lb]);
        fwrite(fileid,[' Target Bodies(1) = 99' lb]);
        fwrite(fileid,[' Name = "Body 1"' lb]);
        fwrite(fileid,[' Equation = 1' lb]);
        fwrite(fileid,[' Material = 1' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);
        
        elseif bodies==2
            
        fwrite(fileid,['Body 1' lb]);
        fwrite(fileid,[' Target Bodies(1) = 56' lb]);
        fwrite(fileid,[' Name = "Body 1"' lb]);
        fwrite(fileid,[' Equation = 1' lb]);
        fwrite(fileid,[' Material = 1' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);
            
        fwrite(fileid,['Body 2' lb]);
        fwrite(fileid,[' Target Bodies(1) = 78' lb]);
        fwrite(fileid,[' Name = "Body 2"' lb]);
        fwrite(fileid,[' Equation = 1' lb]);
        fwrite(fileid,[' Material = 2' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);
                
            
        
        
        
        end
        


        %solve for temp ---------------------------------------------------------------------
        fwrite(fileid,['Solver 1' lb]);
        fwrite(fileid,[' Equation = Heat Equation' lb]);
        fwrite(fileid,[' Procedure = "HeatSolve" "HeatSolver"' lb]);
        fwrite(fileid,[' Variable = Temperature' lb]);
        fwrite(fileid,[' Exec Solver = Always' lb]);
        fwrite(fileid,[' Stabilize = True' lb]);
        fwrite(fileid,[' Bubbles = False' lb]);
        fwrite(fileid,[' Lumped Mass Matrix = False' lb]);
        fwrite(fileid,[' Optimize Bandwidth = True' lb]);
        fwrite(fileid,[' Steady State Convergence Tolerance = 1.0e-5' lb]);
        fwrite(fileid,[' Nonlinear System Convergence Tolerance = 1.0e-7' lb]);
        fwrite(fileid,[' Nonlinear System Max Iterations = 20' lb]);
        fwrite(fileid,[' Nonlinear System Newton After Iterations = 3' lb]);
        fwrite(fileid,[' Nonlinear System Newton After Tolerance = 1.0e-3' lb]);
        fwrite(fileid,[' Nonlinear System Relaxation Factor = 1' lb]);
        fwrite(fileid,[' Linear System Solver = Iterative' lb]);
        fwrite(fileid,[' Linear System Iterative Method = BiCGStab' lb]);
        fwrite(fileid,[' Linear System Max Iterations = 500' lb]);
        fwrite(fileid,[' Linear System Convergence Tolerance = 1.0e-10' lb]);
        fwrite(fileid,[' BiCGstabl polynomial degree = 2' lb]);
        fwrite(fileid,[' Linear System Preconditioning = Diagonal' lb]);
        fwrite(fileid,[' Linear System ILUT Tolerance = 1.0e-3' lb]);
        fwrite(fileid,[' Linear System Abort Not Converged = False' lb]);
        fwrite(fileid,[' Linear System Residual Output = 1' lb]);
        fwrite(fileid,[' Linear System Precondition Recompute = 1' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);




        %save some output files (for paraview ect.)------------------------------------
        fwrite(fileid,['Solver 3' lb]);
        fwrite(fileid,[' Equation = Result Output' lb]);
        fwrite(fileid,[' Procedure = "ResultOutputSolve" "ResultOutputSolver"' lb]);
        fwrite(fileid,[' Output File Name =',name lb]);
        fwrite(fileid,[' Output Format = Vtk' lb]);
        fwrite(fileid,[' Binary Output = False' lb]);
        fwrite(fileid,[' Exec Solver = After all' lb]);  
        fwrite(fileid,[' Scalar Field 1 = Temperature' lb]);
        fwrite(fileid,[' Vector Field 2 = Temperature Grad' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);

        
       
  
        
        
        
  
        fwrite(fileid,['Solver 2' lb]);
        fwrite(fileid,[' Equation = Flux and Gradient' lb]);
        fwrite(fileid,[' Target Variable = Temperature' lb]);
        fwrite(fileid,[' Flux Coefficient = String "Heat Conductivity"' lb]);
        fwrite(fileid,[' Procedure = "FluxSolver" "FluxSolver"' lb]);
        fwrite(fileid,[' Calculate Grad = True' lb]);
       % fwrite(fileid,[' Calculate Flux = True' lb]);
        fwrite(fileid,[' Exec Solver = Always ' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);
        
       

        
        %{ 

        %apply some operators on the data and save results---------------------
        % max, min, max abs, min abs, mean, variance, deviation
        % int mean, int variance
        % volume
        % boundary sum, boundary dofs, boundary mean, boundary max, boundary min,
        % boundary max abs, boundary min abs, area, boundary int, boundary int mean.
        % diffusive energy, convective energy, potential energy
        % diffusive flux, convective flux, boundary int, boundary int mean, area
        % ...
        fwrite(fileid,['Solver 2' lb]);
        fwrite(fileid,['Equation = SaveScalars' lb]);
        fwrite(fileid,['Exec Solver = After All' lb]);
        fwrite(fileid,['Procedure = "SaveData" "SaveScalars"' lb]);
        %fwrite(fileid,['Save Flux = True' lb]);
        fwrite(fileid,['Filename = "'  project_path filesep 'operator_results.dat"' lb]);
        fwrite(fileid,['Variable 1 = potential loads' lb]);
        fwrite(fileid,['Operator 1 = boundary sum' lb]);
        %fwrite(fileid,['Variable 2 =coordinate 1' lb]);
        %fwrite(fileid,['Save Coordinates(1,2) = 0 0' lb])
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);

        %}

        %{
        %save some line data--------------------------------------------------------------------
        fwrite(fileid,['Solver 4' lb]);
        fwrite(fileid,['Exec Solver = After Simulation' lb]);
        fwrite(fileid,['Equation = SaveLine' lb]);
        fwrite(fileid,['Procedure = "SaveData" "SaveLine"' lb]);
        fwrite(fileid,['Filename = "'  project_path filesep 'line.dat"' lb]);
        fwrite(fileid,['Polyline coordinates(2,2) = 0 150e-6 6200e-6 150e-6' lb]);
        %fwrite(fileid,['Save Flux = True' lb]);
        %fwrite(fileid,['Variable 1 = Potential' lb]);
        %fwrite(fileid,['Variable 2 = coordinate 1' lb]);
        %fwrite(fileid,['Variable 3 = coordinate 2' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);

        %}

        %equation-------------------------------------------------------------------------------------
        fwrite(fileid,['Equation 1' lb]);
        fwrite(fileid,[' Name = "Equation 1"' lb]);
        fwrite(fileid,[' Active Solvers(3) = 1 2 3 ' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);
       
        
        
        
        
        
        
        
       

        %material properties-----------------------------------------------------------------------
        fwrite(fileid,['Material 1' lb]);
        fwrite(fileid,[' Name = "Material 1"' lb]);
        fwrite(fileid,[' Heat Conductivity = ',mat2str(mat1) lb]);
        fwrite(fileid,[' Density = 1' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);
        
        fwrite(fileid,['Material 2' lb]);
        fwrite(fileid,[' Name = "Material 2"' lb]);
        fwrite(fileid,[' Heat Conductivity = ',mat2str(mat2) lb]);
        fwrite(fileid,[' Density = 1' lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid, lb);
        
        
        
        
        
         
        
        
        %boundary conditons------------------------------------------------------------------
if nboundaries==1
    
    
        %Condition 1
        fwrite(fileid,['Boundary Condition 1' lb]);
        fwrite(fileid,['  Name = "BoundaryCondition 1"' lb]);
        fwrite(fileid,['  Target Boundaries(',mat2str(length(boundary1)),') = ',num2str(boundary1) lb]);
        fwrite(fileid,['  Temperature = ',mat2str(temp1) lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid,lb);
        
elseif nboundaries==2
    
    %Condition 1
        fwrite(fileid,['Boundary Condition 1' lb]);
        fwrite(fileid,['  Name = "BoundaryCondition 1"' lb]);
        fwrite(fileid,['  Target Boundaries(',mat2str(length(boundary1)),') = ',num2str(boundary1) lb]);
        fwrite(fileid,['  Temperature = ',mat2str(temp1) lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid,lb);
    
  
        %Condition 2
        fwrite(fileid,['Boundary Condition 2' lb]);
        fwrite(fileid,['  Name = "BoundaryCondition 2"' lb]);
        fwrite(fileid,['  Target Boundaries(',mat2str(length(boundary2)),') = ',num2str(boundary2) lb]);
        fwrite(fileid,['  Temperature = ',mat2str(temp2) lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid,lb);
        
elseif nboundaries==3
    
    %Condition 1
        fwrite(fileid,['Boundary Condition 1' lb]);
        fwrite(fileid,['  Name = "BoundaryCondition 1"' lb]);
        fwrite(fileid,['  Target Boundaries(',mat2str(length(boundary1)),') = ',num2str(boundary1) lb]);
        fwrite(fileid,['  Temperature = ',mat2str(temp1) lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid,lb);
    
  
        %Condition 2
        fwrite(fileid,['Boundary Condition 2' lb]);
        fwrite(fileid,['  Name = "BoundaryCondition 2"' lb]);
        fwrite(fileid,['  Target Boundaries(',mat2str(length(boundary2)),') = ',num2str(boundary2) lb]);
        fwrite(fileid,['  Temperature = ',mat2str(temp2) lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid,lb);
        
        
        %Condition 3
        fwrite(fileid,['Boundary Condition 3' lb]);
        fwrite(fileid,['  Name = "BoundaryCondition 3"' lb]);
        fwrite(fileid,['  Target Boundaries(',mat2str(length(boundary3)),') = ',num2str(boundary3) lb]);
        fwrite(fileid,['  Temperature = ',mat2str(temp3) lb]);
        fwrite(fileid,['End' lb]);
        fwrite(fileid,lb);
        
       

        %close file
        fclose(fileid);

        %generate start info for elmersolver.exe=================================================================================
        %define path
        solverinfo_path =  [project_path filesep 'ELMERSOLVER_STARTINFO'];

        %delete old file if it exists
        if isequal(exist(solverinfo_path,'file'),2)
          try
           delete(solverinfo_path);
          catch
          end
        end

           %open file
        fileid = fopen(solverinfo_path,'w');

        %write file
        fwrite(fileid,['case.sif' lb]);
        fwrite(fileid,['1' lb])
        fwrite(fileid, lb)

        %close file
        fclose(fileid);


       

       

end