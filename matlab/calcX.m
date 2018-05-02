function [h] = calcX(I, c, s, x, y, f)
%CALCX 
%   

%% 

n_x = length(x);
n_y = length(y);
err_min = inf;
h_min = zeros(2,n_x);
h_i = zeros(2,n_x);
h_i(1,:) = x;
for i = 0:n_y^n_x-1
    b = dec2base(i, 2, n_x);
    for j = 1:n_x
        h_i(2,j) = y(str2num(b(j)) + 1);
    end
    E = calcE(c, s, h_i, f);
    err = norm(abs(E)-I);
    if err < err_min
        err_min = err;
        h_min = h_i;
    end
end
h = h_min

end
