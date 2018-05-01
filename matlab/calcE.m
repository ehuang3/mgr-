function [E] = calcE(c, s, x, f)
%CALCE 
%

%% 

[~,n_c] = size(c);
[~,n_x] = size(x);
E = complex(zeros(n_c,1));
for i = 1:n_c
    for j = 1:n_x
        dist = norm(x(:,j)-s) + norm(c(:,i)-x(:,j));
        E(i) = E(i) + exp(1i*2*pi/f*dist);
    end
end

end