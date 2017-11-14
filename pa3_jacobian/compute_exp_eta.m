% compute transformation matrix
function exp_eta = compute_exp_eta(w, theta, q)
exp_w = compute_exp_w(w, theta);
p = (eye(3) - exp_w) * q;
exp_eta = [exp_w, p;
           0, 0, 0, 1];
end