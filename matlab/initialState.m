function [ y0 ] = initialState( model, vm )

if (size(model.rs, 2) == 3)
  vars = [1; vm/100; tanh((vm+20)/50.0)];
else
  vars = [1, vm/100];
end

y0 = exp(model.rs * vars);
y0 = y0 / sum(y0);


end
