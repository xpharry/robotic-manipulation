% compute the rotation matrix by w and theta
function exp_w = compute_exp_w(w, theta)
w_hat = skew(w);
exp_w = eye(3) + w_hat * sin(theta) + w_hat * w_hat * (1 - cos(theta));
end
