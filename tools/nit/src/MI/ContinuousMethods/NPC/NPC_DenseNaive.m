function Ker_grid = NPC_DenseNaive(B,data,grid)

    M=5; %%% split the data for memory issues
    
    %%%% for the mex file:
    S=1:size(grid,2);
    Sm=size(grid,2)/M;
    
    if rem(size(S,2),M)~=0
        error('choose knots a product of 5  "or change DneseNaive code"')
    end
    
    if size(data,1)==1
        Ker_gridS=zeros(3,size(grid,2));
    end
    if size(data,1)==2
        Ker_gridS=zeros(5,size(grid,2));
    end
    
    for i=1:M
        SS=(i-1)*Sm+1:i*Sm;
        if i==1            
            
            if size(data,1)==2
                if ispc
                    Ker_gridS = NPC_DenseNaiveMatwin(B,data,grid(:,SS));
                else
                    Ker_gridS = NPC_DenseNaiveMat(B,data,grid(:,SS));
                end
            end
            
        else
            if size(data,1)==2
                
                if ispc
                    Ker_gridS = cat(2,Ker_gridS,NPC_DenseNaiveMatwin(B,data,grid(:,SS)));
                else
                    Ker_gridS = cat(2,Ker_gridS,NPC_DenseNaiveMat(B,data,grid(:,SS)));
                end
                
            end
        end
    end

    if size(data,1)==2
        Ker_grid = mat2cell(Ker_gridS,ones(1,5),size(grid,2))';
    end
    
    
    

