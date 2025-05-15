function [A, B, C, D, F, G] = generate_random_lti_faulty_system(nx,nu,ny,nf,nz)
%GENERATE_RANDOM_LTI_FAULTY_SYSTEM Creates a stable LTI system with given
%dimensions of state, input, output, fault, and fault-to-output zero
%spaces.
%   The system is in discrete time, with dt=1. (A,C) is observable.
%   (A,B) is a.s. controllable, but no control is done for this.
%   This code assumes nf < ny

% Generate A
pole_seed = 2*rand(nx,1)-1;
poles = [];
i = 1;
while i <= nx
    if rand() > 0.5  || i == nx
        poles(i) = pole_seed(i);
    else  % complex pair
        poles(i) = pole_seed(i)*(exp(1i*pi*pole_seed(i+1)));
        poles(i+1) = poles(i)';
        i = i+1;
    end
    i = i+1;
end

den = poly(poles);
A = rot90(compan(den),2);

% B and D matrices are simply random
B = randn(nx,nu);
D = randn(ny,nu);

% C is of the form [Q 0], where Q is orthogonal
[q, ~] = qr(randn(ny,ny));
C = zeros(ny,nx);
C(1:ny, 1:ny) = q;

% Fun fact, this choice gives a nice structure to the singular values of
% the observability matrix: 1s, sqrt(2)s, sqrt(3)s...

% For F and G, we have to place zeros as requested. This is a bit trickier.
% This code considers relatively simple cases for small nf and nz
if nz > 2*nf
    error('Code does not support more than twice more zeros than faults')
end

z = randn(nz,1);  % For now only real zeros
S = {};
for i = 1:length(z)
    S{i} = orth([A - z(i)*eye(nx); C]);
end

F = [];
G = [];

idx_f = 1;
idx_z = 1;
rem_f = nf;
rem_z = nz;

if nf < nz  % Need to create a column in the intersection
    while true
        rem_f = rem_f - 1;
        rem_z = rem_z - 2;
        Sf = common_range(S{idx_z:idx_z+1});
        ns = size(Sf,2);
        v = randn(ns,1);
        F(:,idx_f) = Sf(1:nx,:)*v;
        G(:,idx_f) = Sf(nx+1:end,:)*v;
        if length(tzero(ss(A,F,C,G),1e-10)) == nz
            % Sometimes this part introduces more zeros than planned
            F(:,idx_f+1:nf) = randn(nx,rem_f);
            G(:,idx_f+1:nf) = randn(ny,rem_f);
            return;
            % TODO: it would be nice to understand why and have more
            % control of this process. Beware: one can lose his/her life
            % trying to master MIMO zeros.
        end
        idx_f = idx_f + 1;
        idx_z = idx_z + 2;
        if rem_z <= rem_f
            break;  % debug... occasionally more zeros are introduced
        end
    end
end

% Now we can create one column per column space
for i = idx_z:length(z)
    Si = S{i};
    ns = size(Si,2);
    v = randn(ns,1);
    F(:,idx_f) = Si(1:nx,:)*v;
    G(:,idx_f) = Si(nx+1:end,:)*v;
    idx_f = idx_f + 1;
    rem_f = rem_f - 1;
end

% Finally, fill in the remaining columns randomly
F(:,idx_f:nf) = randn(nx,rem_f);
G(:,idx_f:nf) = randn(ny,rem_f);

    










