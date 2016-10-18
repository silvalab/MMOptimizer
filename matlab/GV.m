function [ vm, gv ] = GV( model )

vm = -120:10:20;
v0 = -120.0;

stepsize = 0.01;
duration = 2.0;
nSteps = ceil(duration / stepsize);

G = zeros(nSteps, length(vm));
y0 = initialState(model, v0);
G(1, :) = model.C * y0;

for i = 1:length(vm)
    A = transitionMatrix(model, vm(i));
    E = expm(A * stepsize); y = y0;
    for j = 2:nSteps
        y = E * y; G(j, i) = model.C * y;
    end
end

gv = max(G); gv = gv / max(gv);

end


