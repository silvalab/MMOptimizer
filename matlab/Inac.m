function [ vm, gv ] = Inac( model )

vm = -120:10:40;
v0 = -120.0;
v1 = 20;
dt = 200;

stepsize = 0.01;
duration = 2.0;
nSteps = ceil(duration / stepsize);

G = zeros(nSteps, length(vm));
y0 = repmat(initialState(model, v0)', length(vm), 1);


for i = 1:length(vm)
    A = transitionMatrix(model, vm(i));
    E = expm(A * dt);
    y0(i, :) = E * y0(i, :)';
end

G(1, :) = y0 * model.C';
A = transitionMatrix(model, v1);
E = expm(A * stepsize);

for i = 1:length(vm)
    y = y0(i, :)';
    for j = 2:nSteps
        y = E * y; G(j, i) = model.C * y;
    end
end

gv = max(G); gv = gv / max(gv);
    

end


