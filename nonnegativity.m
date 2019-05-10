function [out]=nonnegativity(in, options)
% this is a very simple and effective nonnegativity constraint for
% implementation into a variety of programs.
%Options - 1 - steps through the variables given and replaces negative
%values with 0s.
%Options - 2 - steps through the variable given and if all the values for a
%column are negative, multiplies it by a negative 1.
if options == 1
    [M N]=size(in);
    for i = 1:(N)
            for j = 1:M
                if in(j,i) < 0;
                    in(j,i) = 0;
                end
            end
    end
end
count = 0;
if options == 2
    [M N]=size(in);
    for i = 1:(N)
        for j = 1:M
            if in(j,i)<0;
                count=count+1;
            end
        end
        if count==M
            in(:,i)=in(:,i)*-1;
            count=0;
        end
    end
end
if options == 3
    [M N]=size(in);
    for i = 1:(N)
            for j = 1:M
                if in(j,i) < 0;
                    in(j,i) = in(j,i)*-1;
                end
            end
    end
end
[out]=in;