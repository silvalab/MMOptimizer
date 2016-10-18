
% replace the values of graph, rs, rk, G to those found in the .model file

graph = [ ...
    2, 3
    0, 2
    1, 0
    3, 1
    2, 4
    3, 5
    5, 6 ];

rs = [ ...
 -0.1108         -0.9327          1.0234
  5.0274          5.3867         -3.7586
 -3.0856         -6.6859          0.9247
 -0.0935          0.0499         -0.1692
 -5.6538         -4.6885         -6.6175
  1.4184          0.3089          4.9205
  0.3847         -0.6370         -1.1798 ];

rk = [ ...
  -3.3356         -3.3862          2.7316
  1.0967         -0.0264         -1.7839
 -2.2824          4.1722          0.1143
 -3.2755         -1.7381         -1.8883
 -1.5453         -7.2068          2.3778
 -1.8720          2.1327         -0.7924
 -1.0977         -1.8195          1.0538 ];

G = [ 1       0       0       0       0       0       0 ];



ic = incidenceMatrix(graph);
model = struct('ic', ic, 'rs', rs, 'rk', rk, 'C', G);

figure;
hold on;

[vm, gv] = GV(model); 
plot(vm, gv);

[vm, inac] = Inac(model);
plot(vm, inac);

figure; 
hold on;

[t, y] = kinetics(model);
plot(t, y);



