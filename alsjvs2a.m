function [res, Stats, purenbest, crossSecbest, bestiter, Xo1] = alsjvs2a(datacube, spectra, X, Y, C);
%Ok, here is how it works, datacube is the original spectra for the region.
%It is converted into a 2D matrix, spectra is the pure concentration
%profiles that you have already calculated some how. X is the number of
%wavelengths in the datacube, or its third dimension. Y is the number of
%pure components in spectra.
%disp ('% ***********************************************')
%disp ('% What constraints would you like to use in this*')
%disp ('% session of the program?                       *')
%disp ('% 1.Unimodial with non-negativity               *')
%disp ('% 2.Unimodial/non-negativity with Sigmoidal peak*')
%disp ('%   detection.                                  *')
%disp ('% 3.Shape constraints limiting the peaks to a   *')
%disp ('%   combination of gaussian and lorrentzian     *')
%disp ('%   peaks.                                      *')
%disp ('% 4.Both options 1 and 3                        *')
%disp ('% 5.Both options 2 and 3                        *')
%disp ('% 6.Shape constraints limiting the peaks to a   *')
%disp ('%   combination of gaussian and lorrentzian     *')
%disp ('%   peaks with sigmoidal peak detection.        *')
%disp ('%   There are no hard width constraints.        *')
%disp ('% 7.Shape constraints limiting the peaks to a   *')
%disp ('%   combination of gaussian and lorrentzian     *')
%disp ('%   peaks with sigmoidal peak detection.        *')
%disp ('%   Hard width constraints.                     *')
%disp ('% 8.Shape constraints limiting the peaks to a   *')
%disp ('%   combination of gaussian and lorrentzian     *')
%disp ('%   peaks with sigmoidal peak detection.        *')
%disp ('%   Hard width constraints. Last component is   *')
%disp ('% 	made up of random numbers and added to the  *')
%disp ('%   data matrix, it is unconstrained.           *')
%disp ('% 9.Shape constraints limiting the peaks to a   *')
%disp ('%   combination of gaussian and lorrentzian     *')
%disp ('%   peaks with sigmoidal peak detection.        *')
%disp ('%   Hard width constraints. Last component is   *')
%disp ('% 	made up of random numbers and user input,   *')
%disp ('%   it is unconstrained.                        *')
%disp ('% Please make your selection and press enter.   *')
%disp ('% ***********************************************')
%C = input('Option: ');

for i = 1:X;
    [data(:,i)]=datacube(:,2,i);
    [dataNorm(:,i)]=((data(:,i))/(norm(data(:,i))));
    [dataN(i,:)]=(data(:,i)');
end

[dataN1]=[data];
[dataN2]=[data];
%The data should be in a 3 dimensional cube, here we make a 2D array.We
%then normalize it.

for i = 2:Y
    [pure(:,(i-1))]=spectra(:,i);
    [pureN(:,(i-1))]=pure(:,(i-1))/norm(pure(:,(i-1)));
end
if C==8;
    [pureN(:,Y)]=randn(size(pureN(:,1)))*1000;
    Y=Y+1;
end
[pureN1]=[pureN];
[pureN2]=[pureN];
%the spectra matrix, your pure components, should have the first column
%being shift data. As such we do not include this in the fit. This also
%normalizes the spectra as we go.

res=2000000000000000;
resbest=2000000000000000000000000000;
[counts]=1;
[bestiter]=0;
crossSec2=zeros(X,(Y-1));
crossSec=zeros(X,(Y-1));
[Xo1]=zeros((Y-1),4);
[Xo2]=zeros((Y-1),4);
[Xpos]=zeros((Y-1),1);
[Stats]=ones(3,1);

while (counts < 1000);
    if  ((Stats(1,counts)<(10^-6)) && (Stats(2,counts)<(10^-6)) && (Stats(3,counts)<(10^-6)))
        break;
    end;
        dataN1=dataN2;
        crossSec2=crossSec;
        pureN2=pureN1;
        [crossSec]=((dataN(:,:)*pureN1(:,:))/(pureN1(:,:)'*pureN1(:,:)));
        for i = 1:(Y-1)
            for j = 1:X
                if crossSec(j,i) < 0;
                    crossSec(j,i) = 0;
                end
            end
        end
        [pureN1]=(dataN(:,:)'*crossSec(:,:))/(crossSec(:,:)'*crossSec(:,:));
        for i=1:(Y-1);
        [I]=find(pureN1(:,i)<0);
        temp=isempty(I);
        if length(I) == length(pureN1(:,i));
            pureN1(:,i)=pureN1(:,i)*-1;
        elseif temp==0        
            for k=1:length(I)
                pureN1(I(k),i)=0;
            end
        end
        end
        if C~=8 || C~=9;
           [Xo2, I, J, Loc]=sigpeak(pureN1, (Y-1), counts, Xo2);
        end
         if C==8 || C==9;
           [Xo2, I, J, Loc]=sigpeak(pureN1, (Y-2), counts, Xo2);
        end
           for i = 1:(Y-1);
            [purErr]=dataN-crossSec*pureN1'+crossSec(:,i)*pureN1(:,i)';
            [purCom]=[(crossSec(:,i)'*purErr)*[crossSec(:,i)'*crossSec(:,i)]^-1]';
            [pureComS(:,i)]=[purCom'];
            if C==1 || C==4 ;
               [pureN1(:,i)]=junimod(purCom',1,C,Loc(i,1));
            end

            if C==2 || C==5;
               [pureN1(:,i)]=junimod(purCom',1,C,Loc(i,1));
            end

            if C==3 ;
               [pureN1(:,i),Xo1(i,:)]=shapecon(pureN1(:,i),(1),(Xo1(i,:)),counts,(purCom));
            end
            
            if C==4 || C==5 ;
                [pureN1(:,i),Xo1(i,:)]=shapecon(pureN1(:,i),(1),(Xo1(i,:)),counts,(purCom));
            end

            if C==6 ;
                 [pureN1(:,i),Xo1(i,:)]=sigmoidshape(pureN1(:,i),purCom,(1),Xo1(i,:),counts,Xo2(i,:),I);
               
            end
            if C==7 || C==8;
                if C==8 && i==(Y-1)
                    break
                end
                [pureN1(:,i),Xo1(i,:)]=sigmoidshapehard(pureN1(:,i),purCom,(1),Xo1(i,:),counts,Xo2(i,:),I);            
            end
             if C==9
                if C==9 && i==(Y-1)
                    break
                end
                [pureN1(:,i),Xo1(i,:)]=sigmoidshapehard(pureN1(:,i),purCom,(1),Xo1(i,:),counts,Xo2(i,:),I);            
            end
        end
        [dataN2]=pureN1*crossSec';
        [err]=sum(abs(dataN2-dataN1)).^2;
        [res]=sum(abs(err));
            [resbest]=[res];
            [bestiter]=counts;
            [purenbest]=[pureN1];
            [crossSecbest]= [crossSec];
            counts = counts + 1;
            [count]=counts;
            
    [Stats(1,counts)]=[(pureN1(:)/norm(pureN1(:)))-(pureN2(:)/norm(pureN2(:)))]'*[pureN1(:)/norm(pureN1(:))-pureN2(:)/norm(pureN2(:))];
    [Stats(2,counts)]=[crossSec(:)/norm(crossSec)-crossSec2(:)/norm(crossSec2)]'*[crossSec(:)/norm(crossSec)-crossSec2(:)/norm(crossSec2)];
    [Stats(3,counts)]=[dataN2(:)/norm(dataN2(:))-dataN1(:)/norm(dataN1(:))]'*[dataN2(:)/norm(dataN2(:))-dataN1(:)/norm(dataN1(:))];
end