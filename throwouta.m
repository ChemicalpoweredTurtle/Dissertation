function [Output_B] = throwouta(Matrix, Low, H, X);
        Z=1;  
    while Z<(X+1)
    Matrix_Shift(:,Z)=Matrix(:,1,Z);
    Matrix_Data(:,Z)=Matrix(:,2,Z);    
    Z=Z+1;
    end;
    for i=1:X
        [J]=find(Matrix_Shift(:,i)>Low); 
        [K]=find(Matrix_Shift(:,i)>H);
        A=Matrix_Shift(J(1):K(1),i);
        [L]=find(Matrix_Shift(:,1)>Low);
        [M]=find(Matrix_Shift(:,1)>H);
        B=Matrix_Shift(L(1):M(1),1);
        C=Matrix_Data(J(1):K(1),i);
        while length(C)>length(B) 
            temp=B-A(1:length(B));
            [I]=find(temp>0.5);
            if isempty(I)
                A=[A(1:length(A)-1)];
                C=[C(1:length(C)-1)];
            else
                A=[A(1:I(1)-1);A((I(1)+1):length(A))];
                C=[C(1:I(1)-1);C((I(1)+1):length(C))];
            end
            if length(A)==length(B)
                break;      
            end;
        end;
    [IntensityData(:,i)]=C;
    [ShiftData(:,i)]=A; 
    Output{i}=[ShiftData(:,i),IntensityData(:,i)];
    end;
    for i=1:X
        Output_B(:,:,i)=Output{i};
    end
