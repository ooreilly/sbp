p=4;
n=18;
if length(argv()) > 0
        ops = sprintf('resources/%s/', argv(){1});
else
        ops = 'resources/build3/';
end


[Dp,Dm,Hip,Him,Pp,Pm,xp,xm, bp] = sbp_operators(ops, '42', n);

h = xp(2) - xp(1);

Hp = inv(Hip);
Hm = inv(Him);

Qp = Hp*Dp;
Qm = Hm*Dm;
Bp = Qp*0;
Bp(1,1) = -1;
Bp(end,end) = 1;


% SBP-property
tol = 1e-12;
assert(max(max(abs(Qp + Qm' - Bp))) < tol);
assert(max(max(abs(Hp*Pp - (Hm*Pm)'))) < tol);

fprintf('Norm of interpolation operators\n');
fprintf('||PpPm||_2 =  %g \n',norm(Pp*Pm));
fprintf('||PmPp||_2 =  %g \n',norm(Pm*Pp));


zp = zeros(size(Dp,1), size(Dm,2));
zm = zeros(size(Dm,1), size(Dp,2));
A = [zp, Dp - Hip*Bp; 
     Dm, zm];
fprintf('\n');
fprintf('Spectral radius: %g \n',max(abs(eig(A*h))));

% Accuracy
fprintf('\n');
fprintf('Error in differentiation\n');
for i=0:p
    err = Dp*xm.^i - i*xp.^max((i-1),0);
    fprintf('p = %d   Dp max err: %3g \n', i, (max(abs(err(1:bp(1))))));
    err = Dm*xp.^i - i*xm.^max((i-1),0);
    fprintf('p = %d   Dm max err: %3g \n', i, max(abs(err(1:bp(2)))));
end
fprintf('\n');
xm = xm - 0.5;
xp = xp - 0.5;
fprintf('Error in interpolation\n');
for i=0:p
    err = Pp*xm.^i - xp.^i;
    fprintf('p = %d   Pp max err: %3g \n', i, max(abs(err(1:bp(2)))));
    err = Pm*xp.^i - xm.^i;
    fprintf('p = %d   Pm max err: %3g \n', i, max(abs(err(1:bp(1)))));
end

fprintf('\n');
fprintf('Error in interpolation and differentiation\n');
Dcp1 = Dp*Pm;
Dcp2 = Pp*Dm;
Dcm1 = Dm*Pp;
Dcm2 = Pm*Dp;
for i=0:p
    err = Dcp1*xp.^i - i*xp.^max((i-1),0);
    fprintf('p = %d   Dp*Pm max err: %3g \n', i, (max(abs(err(1:bp(1))))));
end
fprintf('\n');
for i=0:p
    err = Dcp2*xp.^i - i*xp.^max((i-1),0);
    fprintf('p = %d   Pp*Dm max err: %3g \n', i, (max(abs(err(1:bp(1))))));
end
fprintf('\n');
for i=0:p
    err = Dcm1*xm.^i - i*xm.^max((i-1),0);
    fprintf('p = %d   Dm*Pp max err: %3g \n', i, max(abs(err(1:bp(2)))));
end
fprintf('\n');
for i=0:p
    err = Dcm2*xm.^i - i*xm.^max((i-1),0);
    fprintf('p = %d   Pm*Dm max err: %3g \n', i, max(abs(err(1:bp(2)))));
end


