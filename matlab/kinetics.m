function [ t, y ] = kinetics( model )

v0 = -120;
vs = 0;

dt = 0.02;
t1 = 5.0;

y0 = initialState(model, v0);
y = zeros(length(0:dt:t1), length(y0));
y(1, :) = y0;

A = transitionMatrix(model, vs);
E = expm(A * dt);

nSteps = length(dt:dt:t1);
for i = 1:nSteps
    y(i+1, :) = (E * y(i, :)')';
end

t = 0:dt:t1;

end


