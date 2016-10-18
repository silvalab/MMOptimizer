function [ ic ] = incidenceMatrix( G )

nNodes = max(G(:)) + 1;
nEdges = length(G);

ic = zeros(nNodes, 2*nEdges);
G = G + 1;

for i = 1:nEdges
    ic(G(i, 1), 2*(i-1) + 1) = -1;
    ic(G(i, 2), 2*(i-1) + 1) = 1;
    ic(G(i, 1), 2*i) = 1;
    ic(G(i, 2), 2*i) = -1;
end

end

