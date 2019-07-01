function [Dp, Dm, Hip, Him, Pp, Pm, xp, xm, bp] = sbp_operators(path, order, n,
                                                            use_sparse)
% function [Dp, Dm, Hip, Him, Pp, Pm, xp, xm, bp] = operators(path, order, n, 
%                                                             use_sparse)
% 
% Arguments:
%  path : path to directory containing operators Dxx, Hxx, Pxx, etc
%  order : String that specifies the order. For example, order = '42'
%  n : Number of grid points.
% 
% Returns:  
%  Dp : Derivative matrix of size (n+1) x (n+2)
%  Dm : Derivative matrix of size (n+2) x (n+1)
%  Pp : Interpolation matrix of size (n+1) x (n+2)
%  Pm : Interpolation matrix of size (n+2) x (n+1)
%  Hp : Norm matrix of size (n+1) x (n+1)
%  Hm : Norm matrix of size (n+2) x (n+2)
%  xp : Grid of size n + 1 
%  xm : Grid of size n + 2 
%  bp(1) : Number of modified boundary points for Dp
%  bp(2) : Number of modified boundary points for Dm

        if nargin < 4
                use_sparse = 1;
        end


        % Build derivative matrix of size (n + 1) x (n + 2)
        op = load(strcat([path, 'D', order, '.mat']));
        s = op.interior;
        idx = op.shift;
        Dp = spdiags(repmat(s, [n+2 1]), idx, n+2, n+2);
        DL = op.left;
        b = size(DL);
        bp(1) = b(1);
        Dp(1:b(1),1:b(2)) = DL;
        Dp(end,:) = [];
        Dp(end-b(1)+1:end,end-b(2)+1:end) = -fliplr(flipud(DL));

        % Build derivative matrix of size (n + 2) x (n + 1)
        op = load(strcat([path, 'Dhat', order, '.mat']));
        s = op.interior;
        idx = op.shift;
        Dm = spdiags(repmat(s, [n+2 1]), idx, n+2, n+2);
        DL = op.left;
        b = size(DL);
        bp(2) = b(1);
        Dm(1:b(1),1:b(2)) = DL;
        Dm(:,end) = [];
        Dm(end-b(1)+1:end,end-b(2)+1:end) = -fliplr(flipud(DL));

        % Build norm matrix of size (n + 1) x (n + 1)
        op = load(strcat([path, 'H', order, '.mat']));
        Hp = speye(n+1);
        HL = op.left;
        b = size(HL);
        for i=1:b(1)
                Hp(i,i) = HL(i);
                Hp(end-i+1,end-i+1) = HL(i);
        end

        % Build norm matrix of size (n + 2) x (n + 2)
        op = load(strcat([path, 'Hhat', order, '.mat']));
        Hm = speye(n+2);
        HL = op.left;
        b = size(HL);
        for i=1:b(1)
                Hm(i,i) = HL(i);
                Hm(end-i+1,end-i+1) = HL(i);
        end

        % Build interpolation matrix of size (n + 1) x (n + 2)
        op = load(strcat([path, 'P', order, '.mat']));
        s = op.interior;
        idx = op.shift;
        Pp = spdiags(repmat(s, [n+2 1]), idx, n+2, n+2);
        PL = op.left;
        b = size(PL);
        Pp(1:b(1),1:b(2)) = PL;
        Pp(end,:) = [];
        Pp(end-b(1)+1:end,end-b(2)+1:end) = fliplr(flipud(PL));

        % Build interpolation matrix of size (n + 2) x (n + 1)
        op = load(strcat([path, 'Phat', order, '.mat']));
        s = op.interior;
        idx = op.shift;
        Pm = spdiags(repmat(s, [n+2 1]), idx, n+2, n+2);
        PL = op.left;
        b = size(PL);
        Pm(1:b(1),1:b(2)) = PL;
        Pm(:,end) = [];
        Pm(end-b(1)+1:end,end-b(2)+1:end) = fliplr(flipud(PL));

        % Grids for the unit interval 0 <= x <= 1
        h = 1.0/n;
        xp = h*[0:n]';
        xm = h*[0 1/2+0:n n]';  

        % Scale matrices by grid spacing
        Dp = Dp/h;
        Dm = Dm/h;
        Hip = inv(Hp)/h;
        Him = inv(Hm)/h;

        if ~use_sparse
                Dp = full(Dp);
                Dm = full(Dm);
                Hip = full(Hp);
                Him = full(Hm);
        end
end
