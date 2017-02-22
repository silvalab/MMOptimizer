
modelFile = '../snapshots/iter_100000.model';
eval(fileread(modelFile));

model = struct('ic', ic, 'rs', rs, 'rk', rk, 'C', G, 'args', args);

figure(1);
hold on;

[vm, gv] = GV(model); 
plot(vm, gv);

[vm, inac] = Inac(model);
plot(vm, inac);

figure(2); 
hold on;

[t, y] = kinetics(model);
plot(t, y);



