function [Output_B] = throwina(Matrix, Low, H, X);
    Z=1;
    while Z<(X+1)
    Matrix_Shift(:,Z)=Matrix(:,1,X);
    Matrix_Data(:,Z)=Matrix(:,2,X);    
    Z=Z+1;
    end;
for i=1:X
    [L]=find(Matrix_Shift(:,6)>Low);
    [M]=find(Matrix_Shift(:,6)>H);    
    B=Matrix_Shift(L(1)-1:M(1),6);
    [J]=find(Matrix_Shift(:,i)>Low);
    [K]=find(Matrix_Shift(:,i)>H);
    A=Matrix_Shift(J(1)-1:K(1),i);
    C=Matrix_Data(J(1)-1:K(1),i);
    while length(C)<length(B)  
        temp=A-B(1:length(A));
        for j=1:length(A);
            if temp(j) > 0.5 ; 
                A=[A(1:(j-1));mean(A((j-1):j));A(j:length(A))];
                C=[C(1:(j-1));mean(C((j-1):j));C(j:length(C))];
            if length(A)==length(B)
                break;      
            end;
            elseif  temp(j) < -0.5 ; 
               A=[A(1:(j-1));mean(A(j:(j+1)));A(j:length(A))];
               C=[C(1:(j-1));mean(C(j:(j+1)));C(j:length(C))];
            if length(A)==length(B)
                break;      
            end;
            end;
        end
            if length(A)==length(B)
                break;      
            end;
    end;
    
ShiftData(:,i)=A; 
IntensityData(:,i)=C;
 
Output{i}=[A,C];
end;

for G=1:X
    Output_B(:,:,G)=Output{G};
end
