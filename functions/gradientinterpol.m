function     Gradient=gradientinterpol(Gradient_ori,node_iris,pointsTetra,bar)
             %Field interpolation using baricentric coordinates and GPU

                        DT=Gradient_ori(pointsTetra,:);
                        DT1=reshape(DT',3,node_iris,4);
                        DT2=permute(DT1,[3 1 2]);
                        M1= gpuArray(DT2);
                        bar2=reshape(bar',1,4,node_iris);
                        bar22=repmat(bar2,3,1,1);
                        bar222=permute( bar22,[2 1 3]);
                          M2= gpuArray(bar222);
                          Q2= pagefun(@times, M1, M2);
                          Q2= gather(Q2);
                          Q3=sum(Q2,1);
                        Gradient=reshape(Q3,3,node_iris)';
                        
                        gpuDevice(1)


end