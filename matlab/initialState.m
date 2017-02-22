function [ y0 ] = initialState( model, vm )

args = zeros(size(model.args, 1)+1, 1);
args(1) = 1.0;

for i = 1:size(model.args, 1)
    args(i+1) = tanh((vm + model.args(i, 1))/model.args(i, 2));
end


y0 = exp(model.rs * args);
y0 = y0 / sum(y0);


end
