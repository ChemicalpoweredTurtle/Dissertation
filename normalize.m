function [out]=normalize(in, options)
[m n]=size(in);
if options==1
    for i=1:n
        in(:,i)=in(:,i)/sum(in(:,i));
    end
    [out]=in;
end
if options==2
    for i=1:n
        in(:,i)=in(:,i)./(ones(length(in(:,i)),1)*std(in(:,i)));
    end
[out]=in;
end
