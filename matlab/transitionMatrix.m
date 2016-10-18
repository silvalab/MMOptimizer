function [ Q ] = transitionMatrix( model , vm )

ic = model.ic;
rs = model.rs;
rk = model.rk;

if ( size(rs, 2) == 3 )
  args = [1; vm/100; tanh((vm+20)/50.0)];
else
  args = [1; vm /100];
end

N = size(rs, 1);
E = size(rk, 1);

D = zeros(E, 2*E);
D(:, 1:2:end) = eye(E);
D(:, 2:2:end) = -eye(E);
D = [D;abs(D)];

r_vec = exp(D^-1 * [ic(:, 1:2:end)'*rs; rk] * args);
Q = ic * diag(r_vec(:)) * (-min(ic, 0)');

end
