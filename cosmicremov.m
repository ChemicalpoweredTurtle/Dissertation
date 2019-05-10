function [clean]=cosmicremov(data, X);

[m,n] = size(data);
[store]=data;
[C] = 0;
for i = 1:m
    for j = 1:n
        [Ave] = (sum(data(i,:)) - data(i,j))/(n-1);
       if abs(data(i,j))<(abs(Ave)-abs((X*Ave)))
            store(i,j) = Ave;
            C = C+1;
        end
        if  data(i,j)>(Ave+(X*Ave))
            store(i,j) = Ave;
            C = C+1;
        end
    end
end
disp(C);
[clean] = store;
