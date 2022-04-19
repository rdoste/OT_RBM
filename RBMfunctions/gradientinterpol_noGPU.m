function     Gradient=gradientinterpol_noGPU(Gradient_ori,node_iris,pointsTetra,bar)




                        DT=Gradient_ori(pointsTetra,:);
                        DT1=reshape(DT',3,node_iris,4);
                        DT2=permute(DT1,[3 1 2]);
                     
                        bar2=reshape(bar',1,4,node_iris);
                        bar22=repmat(bar2,3,1,1);
                        bar222=permute( bar22,[2 1 3]);
                         
                      
                        
                   
                        
                QX2=zeros(size(DT2));        
              for indx=1:length(bar)
                     QX2(:,:,indx)= times(DT2(:,:,indx), bar222(:,:,indx));
              end
                QX3=sum(QX2,1);        
                GR=QX3(1,:,:);
                Gradient=reshape(GR,3,node_iris)';



end