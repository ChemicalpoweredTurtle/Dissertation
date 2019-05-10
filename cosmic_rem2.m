function [data_new, replaced] = cosmic_rem2(data, alpha, alpha2);
% Tool to remove cosmic spikes observed in long collection time data
data_new = data; 

[n, m]=size(data);

datam = data - ones(n,1)*mean(data,1);
replaced = [];

C = [data'*data]/(n-1);
for i = 1:n;  
    D2(i) = [data(i,:)-mean(data,1)]*inv(C)*[data(i,:)-mean(data,1)]';
end

scale=[0.01:0.01:3];

[Fpdf1] = fpdf(scale,m,n-m);
for i=1:length(Fpdf1);
    [Fpdf1_area_per(i)] = trapz(Fpdf1(1:i))/trapz(Fpdf1);
end

[L1] = find(Fpdf1_area_per > (1-alpha));
T2_crit1 = ([m*(n-1)*(n+1)]/[n*(n-m)])*scale(L1(1));

%figure,hold
%plot([D2',T2_crit1*ones(n,1)])

% Need to set this up to deal with spikes on the ends of the spectrum.


clear D2_2
[I] = find(D2 > T2_crit1);
for i = 1:length (I)
    if I(i) == 1
        C = (datam(I(i):I(i)+4,:)*datam(I(i):I(i)+4,:)')/(m-1);
        for j = 1:m
            D2_2(i,j) = [datam(I(i):I(i)+4,j)-mean(datam(I(i):I(i)+4,:),2)]'*inv(C)*[datam(I(i):I(i)+4,j)-mean(datam(I(i):I(i)+4,:),2)];
        end
    elseif I(i) == n
        C = (data(I(i)-4:I(i),:)*data(I(i)-4:I(i),:)')/(m-1);
        for j = 1:m
            D2_2(i,j) = [datam(I(i)-4:I(i),j)-mean(datam(I(i)-4:I(i),:),2)]'*inv(C)*[datam(I(i)-4:I(i),j)-mean(datam(I(i)-4:I(i),:),2)];
        end
    else
        C = (datam(I(i)-2:I(i)+2,:)*datam(I(i)-2:I(i)+2,:)')/(m-1);
        for j = 1:m
            D2_2(i,j) = [datam(I(i)-2:I(i)+2,j)-mean(datam(I(i)-2:I(i)+2,:),2)]'*inv(C)*[datam(I(i)-2:I(i)+2,j)-mean(datam(I(i)-2:I(i)+2,:),2)];
        end
    end
end


scale2=[0.1:0.1:10];
[Fpdf2] = fpdf(scale2,5,m-5);
for i=1:length(Fpdf2);
    [Fpdf2_area_per(i)] = trapz(Fpdf2(1:i))/trapz(Fpdf2);
end

[L2] = find(Fpdf2_area_per > (1-alpha2));
T2_crit2 = ([5*(m-1)*(m+1)]/[m*(m-5)])*scale2(L2(1));

%figure,hold
%plot([D2_2',T2_crit2*ones(m,1)])

[L M]= find(D2_2>T2_crit2);

for i = 1:length(L);
    temp = find (I(L)==I(L(i)));
    in1 = [1:m]';
    for j=1:length(temp);
        in1 = find (in1~=M(temp(j))); 
    end
    data_new(I(L(i)),M(i))= mean(data(I(L(i)),in1));
    replaced (i,:) = [I(L(i)),M(i)];
end
    
    