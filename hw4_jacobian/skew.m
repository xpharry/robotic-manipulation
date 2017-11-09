% compute the (x)^, namely the skew symmetric matrix
function x_hat = skew(x)
x_hat = [0 -x(3) x(2);
         x(3) 0 -x(1);
         -x(2) x(1) 0];
end