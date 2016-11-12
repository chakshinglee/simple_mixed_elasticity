% Run script for solving linear elasticity equation (fine scale) with mixed 
% finite element method 

% function [errL2u, errL2s] = fine_problem(nx,ny)
% function s11h = fine_problem(kkk,n)
% clear all;
nx = 101;
ny = nx;

hx = 1/(nx-1);
hy = 1/(ny-1);

% lambda = 1*ones(1,nElements);
% mu = 1/2*ones(1,nElements);

% load coef.mat
% lambda = .22/(1+.22)/(1-2*.22)*coef;
% mu = .5/(1+.22)*coef;

nu = .499;
lambda = nu/(1+nu)/(1-2*nu);
mu = .5/(1+nu);

u = @(x)[4*x(:,1).*(1-x(:,1)).*x(:,2).*(1-x(:,2)) -4*x(:,1).*(1-x(:,1)).*x(:,2).*(1-x(:,2))];
f = @(x)[1*ones(size(x,1),1) -1*ones(size(x,1),1)];
% f = @(x)[(x(:,1)<.5)-(x(:,1)>.5) -(x(:,2)<.5)+(x(:,2)>.5)];
% f = @(x)[(2*mu+lambda).*8.*x(:,2).*(x(:,2)-1)+2*mu.*4*x(:,1).*(x(:,1)-1)-(4*lambda+4*mu).*(1-2*x(:,1)).*(1-2*x(:,2)) ...
%     (4*lambda+4*mu).*(1-2*x(:,1)).*(1-2*x(:,2))-(2*mu+lambda).*8.*x(:,1).*(x(:,1)-1)-2*mu.*4*x(:,2).*(x(:,2)-1)];

lambda1 = lambda;
mu1 = mu;

sigma11 = @(x)(2*mu1+lambda1).*4.*(1-2*x(:,1)).*x(:,2).*(1-x(:,2))-lambda1.*4.*x(:,1).*(1-x(:,1)).*(1-2*x(:,2));
sigma22 = @(x)lambda1.*4.*(1-2*x(:,1)).*x(:,2).*(1-x(:,2))-(2*mu1+lambda1).*4.*x(:,1).*(1-x(:,1)).*(1-2*x(:,2));
sigma12 = @(x)2*mu1.*2.*(x(:,1).*(1-x(:,1)).*(1-2*x(:,2))-(1-2*x(:,1)).*x(:,2).*(1-x(:,2)));

lambda = lambda*ones((nx-1)*(ny-1),1);
mu = mu*ones((nx-1)*(ny-1),1);


%% Mesh and topology
Mesh = TProd_Mesh(0:hx:1,0:hy:1);
nElements = size(Mesh.Elements,1);
nCoordinates = size(Mesh.Coordinates,1);

el_node = sparse([1:nElements 1:nElements 1:nElements 1:nElements],Mesh.Elements(:),[ones(1,nElements) 2*ones(1,nElements) 3*ones(1,nElements) 4*ones(1,nElements)]);
node_node = el_node'*el_node;
node_node = node_node-diag(diag(node_node));
node_node(node_node==3)=0;
node_node(node_node==8)=0;
node_node = tril(node_node);
[In,Jn,~]=find(node_node);
nEdges = numel(In);
Mesh.Vert2Edge = sparse(In,Jn,1:numel(In),nCoordinates,nCoordinates);
Mesh.Vert2Edge = Mesh.Vert2Edge + Mesh.Vert2Edge'; 

clear In Jn node_node

% Mesh = add_Edges(Mesh); nEdges = size(Mesh.Edges,1);
% Loc = get_BdEdges(Mesh);
% Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
% Mesh.BdFlags(Loc) = -1;

% plot_Mesh(Mesh,'tas');
% plot_Mesh(Mesh,'tas');

QuadRule_1D = gauleg(0,1,2);
QuadRule = TProd(QuadRule_1D);

%% [Ib, Jb, B] = assembleDivSigma(Mesh);

a11 = (lambda+2*mu)/4./mu./(mu+lambda);
a12 = -lambda/4./mu./(mu+lambda);
a33 = 1/2./mu;

Ib = zeros(nElements*4+nElements*8,1);
Jb = Ib;
Ab = Ib;

% Edge dof
tmp = ones(4,nElements);
tmp(1,:) = -hy;
tmp(2,:) = hy;
tmp(3,:) = -hx;
tmp(4,:) = hx;
Ab(1:nElements*4) = tmp(:); % (-1/hx)*hx*hy and so on
% EdgeLoc = [diag(Mesh.Vert2Edge(Mesh.Elements(:,4),Mesh.Elements(:,1))) diag(Mesh.Vert2Edge(Mesh.Elements(:,2),Mesh.Elements(:,3))) ...
%     diag(Mesh.Vert2Edge(Mesh.Elements(:,1),Mesh.Elements(:,2))) diag(Mesh.Vert2Edge(Mesh.Elements(:,4),Mesh.Elements(:,3)))]';
EdgeLoc = [Mesh.Vert2Edge(Mesh.Elements(:,4)+(Mesh.Elements(:,1)-1)*nCoordinates) Mesh.Vert2Edge(Mesh.Elements(:,2)+(Mesh.Elements(:,3)-1)*nCoordinates) ...
    Mesh.Vert2Edge(Mesh.Elements(:,1)+(Mesh.Elements(:,2)-1)*nCoordinates) Mesh.Vert2Edge(Mesh.Elements(:,4)+(Mesh.Elements(:,3)-1)*nCoordinates)]';
Ib(1:nElements*4) = nCoordinates + EdgeLoc(:);
tmp = [2*(1:nElements)-1;2*(1:nElements)-1;2*(1:nElements);2*(1:nElements)];
Jb(1:nElements*4) = nEdges + nCoordinates + tmp(:);
EdgeLoc = EdgeLoc';

% Vertice dof
tmp = ones(4,nElements);
tmp(1:2,:) = -1;
Ab(nElements*4+1:nElements*4+nElements*4) = tmp(:)*hx; % same reason as for edges
tmp = ones(4,nElements);
tmp([1 4],:) = -1;
Ab(nElements*4+nElements*4+1:end) = tmp(:)*hy;
tmp = Mesh.Elements';
Ib(nElements*4+1:end) = [tmp(:); tmp(:)];
tmp = [1:nElements;1:nElements;1:nElements;1:nElements];
Jb(nElements*4+1:nElements*4+nElements*4) = nEdges + nCoordinates + 2*tmp(:) -1;
Jb(nElements*4+nElements*4+1:end) = nEdges + nCoordinates + 2*tmp(:);
clear tmp;

%% [Im, Jm, M] = assembleMassSigma(Mesh);

Im = zeros(nElements*32,1);
Jm = Im;
Am = Im;

% in each element, each edge has nonzeros with other edges, and each point
% has nonzeros with other points

%
%    4 ___ 3          ___
%     |   |          | 4 |
%     |___|        1 |___| 2
%    1     2           3

Am(1:nElements*10) = hx*hy*[sum(QuadRule.x(:,1).^2.*QuadRule.w)*a11 sum((1-QuadRule.x(:,1)).^2.*QuadRule.w)*a11 ...
    sum(QuadRule.x(:,2).^2.*QuadRule.w)*a11 sum((1-QuadRule.x(:,2)).^2.*QuadRule.w)*a11 ...
    sum(QuadRule.x(:,1).*(1-QuadRule.x(:,1)).*QuadRule.w)*a11 sum(QuadRule.x(:,1).*QuadRule.x(:,2).*QuadRule.w)*a12 ...
    sum(QuadRule.x(:,1).*(1-QuadRule.x(:,2)).*QuadRule.w)*a12 sum((1-QuadRule.x(:,1)).*QuadRule.x(:,2).*QuadRule.w)*a12 ...
    sum((1-QuadRule.x(:,1)).*(1-QuadRule.x(:,2)).*QuadRule.w)*a12 sum(QuadRule.x(:,2).*(1-QuadRule.x(:,2)).*QuadRule.w)*a11];
Am(nElements*10+1:nElements*16) = Am(nElements*4+1:nElements*10); % symmetry

Im(1:nElements*16) = [EdgeLoc(:,2) EdgeLoc(:,1) EdgeLoc(:,4) EdgeLoc(:,3) EdgeLoc(:,2) EdgeLoc(:,2) EdgeLoc(:,2) EdgeLoc(:,1) EdgeLoc(:,1) EdgeLoc(:,4) ...
    EdgeLoc(:,1) EdgeLoc(:,4) EdgeLoc(:,3) EdgeLoc(:,4) EdgeLoc(:,3) EdgeLoc(:,3)] + nCoordinates;
Jm(1:nElements*16) = [EdgeLoc(:,2) EdgeLoc(:,1) EdgeLoc(:,4) EdgeLoc(:,3) EdgeLoc(:,1) EdgeLoc(:,4) EdgeLoc(:,3) EdgeLoc(:,4) EdgeLoc(:,3) EdgeLoc(:,3) ...
    EdgeLoc(:,2) EdgeLoc(:,2) EdgeLoc(:,2) EdgeLoc(:,1) EdgeLoc(:,1) EdgeLoc(:,4)] + nCoordinates;

Am(nElements*16+1:nElements*26) = 2*hx*hy*[sum((QuadRule.x(:,1)+QuadRule.x(:,2)-.5).^2.*QuadRule.w)*a33 sum((QuadRule.x(:,1)-QuadRule.x(:,2)+.5).^2.*QuadRule.w)*a33 ...
    sum((1.5-QuadRule.x(:,1)-QuadRule.x(:,2)).^2.*QuadRule.w)*a33 sum((.5-QuadRule.x(:,1)+QuadRule.x(:,2)).^2.*QuadRule.w)*a33 ...
    sum((QuadRule.x(:,1)+QuadRule.x(:,2)-.5).*(QuadRule.x(:,1)-QuadRule.x(:,2)+.5).*QuadRule.w)*a33 ...
    sum((QuadRule.x(:,1)+QuadRule.x(:,2)-.5).*(1.5-QuadRule.x(:,1)-QuadRule.x(:,2)).*QuadRule.w)*a33 ...
    sum((QuadRule.x(:,1)+QuadRule.x(:,2)-.5).*(.5-QuadRule.x(:,1)+QuadRule.x(:,2)).*QuadRule.w)*a33 ...
    sum((QuadRule.x(:,1)-QuadRule.x(:,2)+.5).*(1.5-QuadRule.x(:,1)-QuadRule.x(:,2)).*QuadRule.w)*a33 ...
    sum((QuadRule.x(:,1)-QuadRule.x(:,2)+.5).*(.5-QuadRule.x(:,1)+QuadRule.x(:,2)).*QuadRule.w)*a33 ...
    sum((1.5-QuadRule.x(:,1)-QuadRule.x(:,2)).*(.5-QuadRule.x(:,1)+QuadRule.x(:,2)).*QuadRule.w)*a33];
Am(nElements*26+1:nElements*32) = Am(nElements*20+1:nElements*26);   % symmetry

Im(nElements*16+1:nElements*32) = [Mesh.Elements(:,3)' Mesh.Elements(:,2)' Mesh.Elements(:,1)' Mesh.Elements(:,4)' Mesh.Elements(:,3)' Mesh.Elements(:,3)' Mesh.Elements(:,3)' ...
    Mesh.Elements(:,2)' Mesh.Elements(:,2)' Mesh.Elements(:,1)' Mesh.Elements(:,2)' Mesh.Elements(:,1)' Mesh.Elements(:,4)' Mesh.Elements(:,1)' Mesh.Elements(:,4)' Mesh.Elements(:,4)'];
Jm(nElements*16+1:nElements*32) = [Mesh.Elements(:,3)' Mesh.Elements(:,2)' Mesh.Elements(:,1)' Mesh.Elements(:,4)' Mesh.Elements(:,2)' Mesh.Elements(:,1)' Mesh.Elements(:,4)' ...
    Mesh.Elements(:,1)' Mesh.Elements(:,4)' Mesh.Elements(:,4)' Mesh.Elements(:,3)' Mesh.Elements(:,3)' Mesh.Elements(:,3)' Mesh.Elements(:,2)' Mesh.Elements(:,2)' Mesh.Elements(:,1)'];

A = sparse([Ib; Jb; Im],[Jb; Ib; Jm], [Ab; Ab; Am]);

% activeDof = setdiff(2:nEdges+nCoordinates+nElements*2,[nCoordinates+[EdgeLoc(1:ny-1,1); EdgeLoc(end-ny+2:end,2); ...
%     EdgeLoc(1:ny-1:(ny-1)*(nx-1),3); EdgeLoc(ny-1:ny-1:(ny-1)*(nx-1),4)];Mesh.Elements(1:ny-1,4); Mesh.Elements(end-ny+2:end,2); Mesh.Elements(1:(ny-1):(ny-1)*(nx-1),1); ...
%     Mesh.Elements(ny-1:(ny-1):(ny-1)*(nx-1),3); nEdges+nCoordinates+[(ny-2)*(nx-1)+1;(ny-1)*(nx-1);(ny-1)]]);
activeDof = 2:nEdges+nCoordinates+nElements*2;

%% Assemble right hand side

F = zeros(2,nElements);
nPt = numel(QuadRule.w);
ScaledQuadRule = QuadRule.x*[hx 0;0 hy];
for i = 1:nPt
    F = F + QuadRule.w(i)'*f([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)])';
end
F = [zeros(nCoordinates+nEdges,1);F(:)*hx*hy];
% U0 = 0*F; U0(nCoordinates+EdgeLoc(kkk,1))=1;
% U0([Mesh.Elements(1:ny-1,4); Mesh.Elements(end-ny+2:end,2); Mesh.Elements(1:(ny-1):(ny-1)*(nx-1),1); ...
%     Mesh.Elements(ny-1:(ny-1):(ny-1)*(nx-1),3)])=([1-(1:ny-1)/(ny-1) -(0:nx-2)/(nx-1) 1-(0:ny-2)/(ny-1) -(1:nx-1)/(nx-1)]+1)/2;
% F = F-A*U0;

%% Compute error

Uh = zeros(nCoordinates+nEdges+nElements*2,1);

errL2u = 0;
errL2s = 0;
for i = 1:nPt
    errL2u = errL2u + QuadRule.w(i)'*sum(sum((u([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
        - [Uh(end-nElements*2+1:2:end) Uh(end-nElements*2+2:2:end)]).^2));
    errL2s = errL2s + QuadRule.w(i)*sum((sigma11([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
        -Uh(nCoordinates+EdgeLoc(:,1))*(1-QuadRule.x(i,1))-Uh(nCoordinates+EdgeLoc(:,2))*QuadRule.x(i,1)).^2);
    errL2s = errL2s + QuadRule.w(i)*sum((sigma22([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
        -Uh(nCoordinates+EdgeLoc(:,3))*(1-QuadRule.x(i,2))-Uh(nCoordinates+EdgeLoc(:,4))*QuadRule.x(i,2)).^2);
    errL2s = errL2s + 2*QuadRule.w(i)*sum((sigma12([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
        -Uh(Mesh.Elements(:,1))*(1.5-QuadRule.x(i,1)-QuadRule.x(i,2))-Uh(Mesh.Elements(:,2))*(.5+QuadRule.x(i,1)-QuadRule.x(i,2)) ...
        -Uh(Mesh.Elements(:,3))*(QuadRule.x(i,1)+QuadRule.x(i,2)-.5)-Uh(Mesh.Elements(:,4))*(.5-QuadRule.x(i,1)+QuadRule.x(i,2))).^2);
end
normL2s = sqrt(errL2s*hx*hy);

% U = u(Mesh.Coordinates(Mesh.Elements(:,1),:)+ones(nElements,2)*[hx/2 0; 0 hy/2])';
normL2u = sqrt(errL2u*hx*hy);

% Uh = U0;
Uh(activeDof) = A(activeDof,activeDof)\F(activeDof);

% U = zeros(2,nElements);
% for i = 1:nElements
%     U(2*i-1:2*i) = u(Mesh.Coordinates(Mesh.Elements(i,1),:)+[hx/2 hy/2]);
% end

errL2u = 0;
errL2s = 0;
for i = 1:nPt
    errL2u = errL2u + QuadRule.w(i)'*sum(sum((u([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
        - [Uh(end-nElements*2+1:2:end) Uh(end-nElements*2+2:2:end)]).^2));
    errL2s = errL2s + QuadRule.w(i)*sum((sigma11([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
        -Uh(nCoordinates+EdgeLoc(:,1))*(1-QuadRule.x(i,1))-Uh(nCoordinates+EdgeLoc(:,2))*QuadRule.x(i,1)).^2);
    errL2s = errL2s + QuadRule.w(i)*sum((sigma22([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
        -Uh(nCoordinates+EdgeLoc(:,3))*(1-QuadRule.x(i,2))-Uh(nCoordinates+EdgeLoc(:,4))*QuadRule.x(i,2)).^2);
    errL2s = errL2s + 2*QuadRule.w(i)*sum((sigma12([Mesh.Coordinates(Mesh.Elements(:,1),1)+ScaledQuadRule(i,1) Mesh.Coordinates(Mesh.Elements(:,1),2)+ScaledQuadRule(i,2)]) ...
        -Uh(Mesh.Elements(:,1))*(1.5-QuadRule.x(i,1)-QuadRule.x(i,2))-Uh(Mesh.Elements(:,2))*(.5+QuadRule.x(i,1)-QuadRule.x(i,2)) ...
        -Uh(Mesh.Elements(:,3))*(QuadRule.x(i,1)+QuadRule.x(i,2)-.5)-Uh(Mesh.Elements(:,4))*(.5-QuadRule.x(i,1)+QuadRule.x(i,2))).^2);
end
errL2s = sqrt(errL2s*hx*hy)/normL2s;

% U = u(Mesh.Coordinates(Mesh.Elements(:,1),:)+ones(nElements,2)*[hx/2 0; 0 hy/2])';
errL2u = sqrt(errL2u*hx*hy)/normL2u;
fprintf('Relative L2 error of the displacement is %f\n',errL2u);
fprintf('Relative L2 error of the stress is %f\n',errL2s);

% plot_BFE(Uh(end-nElements*2+1:2:end),Mesh)

s11h = [Uh(nCoordinates+EdgeLoc(:,1));Uh(nCoordinates+EdgeLoc(:,2));Uh(nCoordinates+EdgeLoc(:,2));Uh(nCoordinates+EdgeLoc(:,1))];
s22h = [Uh(nCoordinates+EdgeLoc(:,3));Uh(nCoordinates+EdgeLoc(:,3));Uh(nCoordinates+EdgeLoc(:,4));Uh(nCoordinates+EdgeLoc(:,4))];
s12h = [1.5*Uh(Mesh.Elements(:,1))+.5*Uh(Mesh.Elements(:,2))-.5*Uh(Mesh.Elements(:,3))+.5*Uh(Mesh.Elements(:,4)); ...
    .5*Uh(Mesh.Elements(:,1))+1.5*Uh(Mesh.Elements(:,2))+.5*Uh(Mesh.Elements(:,3))-.5*Uh(Mesh.Elements(:,4)); ...
    -.5*Uh(Mesh.Elements(:,1))+.5*Uh(Mesh.Elements(:,2))+1.5*Uh(Mesh.Elements(:,3))+.5*Uh(Mesh.Elements(:,4)); ...
    .5*Uh(Mesh.Elements(:,1))-.5*Uh(Mesh.Elements(:,2))+.5*Uh(Mesh.Elements(:,1))+1.5*Uh(Mesh.Elements(:,4))];
% 
% midpt = sum([Mesh.Coordinates(Mesh.Elements(:,1),1) Mesh.Coordinates(Mesh.Elements(:,2),1) Mesh.Coordinates(Mesh.Elements(:,3),1) ...
%     Mesh.Coordinates(Mesh.Elements(:,4),1)],2);
% midpt = [midpt sum([Mesh.Coordinates(Mesh.Elements(:,1),2) Mesh.Coordinates(Mesh.Elements(:,2),2) Mesh.Coordinates(Mesh.Elements(:,3),2) ...
%     Mesh.Coordinates(Mesh.Elements(:,4),2)],2)]/4;
% % s11 = sigma11(midpt);
% % s22 = sigma22(midpt);
% s12 = sigma12(midpt);

% plot_DGBFE(s12h,Mesh)
