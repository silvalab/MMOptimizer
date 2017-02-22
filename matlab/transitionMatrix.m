function [ Q ] = transitionMatrix( model , vm )

ic = model.ic;
rs = model.rs;
rk = model.rk;

args = zeros(size(model.args, 1)+1, 1);
args(1) = 1.0;

for i = 1:size(model.args, 1)
    args(i+1) = tanh((vm + model.args(i, 1))/model.args(i, 2));
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
