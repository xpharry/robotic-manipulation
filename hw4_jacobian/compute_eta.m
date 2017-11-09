% compute transformation matrix
function eta = compute_eta(w, q)
eta = [-cross(w, q);
       w];
end