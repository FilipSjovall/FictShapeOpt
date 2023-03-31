function mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d)# ok
    
    epsimin = 10^(-7);
    raa0 = 0.00001;
    move = 0.5;
    albefa = 0.1;
    asyinit = 0.5;
    asyincr = 1.2;
    asydecr = 0.7;
    eeen = ones(n,1);
    eeem = ones(m,1);
    zeron = zeros(n,1);

    # Calculation of the asymptotes low and upp :
    if iter < 2.5
      low = xval - asyinit*(xmax-xmin);
      upp = xval + asyinit*(xmax-xmin);
    else
      zzz = (xval-xold1).*(xold1-xold2);
      factor = eeen;
        factor = reshape([ (asyincr*(x > 0) + asydecr*(x<0)) for x in zzz[:] ],size( zzz))
      low = xval - factor.*(xold1 - low);
      upp = xval + factor.*(upp - xold1);
      lowmin = xval - 10*(xmax-xmin);
      lowmax = xval - 0.01*(xmax-xmin);
      uppmin = xval + 0.01*(xmax-xmin);
      uppmax = xval + 10*(xmax-xmin);
      low = max.(low,lowmin);
      low = min.(low,lowmax);
      upp = min.(upp,uppmax);
      upp = max.(upp,uppmin);
    end

    # Calculation of the bounds alfa and beta :

    zzz1 = low + albefa*(xval-low);
    zzz2 = xval - move*(xmax-xmin);
    zzz  = max.(zzz1,zzz2);
    alfa = max.(zzz,xmin);
    zzz1 = upp - albefa*(upp-xval);
    zzz2 = xval + move*(xmax-xmin);
    zzz  = min.(zzz1,zzz2);
    beta = min.(zzz,xmax);

    # Calculations of p0, q0, P, Q and b.

    xmami = xmax-xmin;
    xmamieps = 0.00001*eeen;
    xmami = max.(xmami,xmamieps);
    xmamiinv = eeen./xmami;
    ux1 = upp-xval;
    ux2 = ux1.*ux1;
    xl1 = xval-low;
    xl2 = xl1.*xl1;
    uxinv = eeen./ux1;
    xlinv = eeen./xl1;
    
    p0 = zeron;
    q0 = zeron;
    #p0 = max.(df0dx,0);
    #q0 = max.(-df0dx,0);
    p0 = max.(df0dx,0);
    q0 = max.(-df0dx,0);
    pq0 = 0.001*(p0 + q0) + raa0*xmamiinv;
    p0 = p0 + pq0;
    q0 = q0 + pq0;
    p0 = p0.*ux2;
    q0 = q0.*xl2;


    P = max.(dfdx,0);
    Q = max.(-dfdx,0);
    P = sparse(P);
    Q = sparse(Q);
    PQ = 0.001*(P + Q) + raa0*eeem*xmamiinv';
    P = P + PQ;
    Q = Q + PQ;
    P = P * sparse(Diagonal(ux2[:,1]))
    Q = Q * sparse(Diagonal(xl2[:,1]))#spdiags(xl2,0,n,n);
    b = P*uxinv + Q*xlinv .- fval ;
    #Solving the subproblem by a primal-dual Newton method
    xmma,ymma,zmma,lam,xsi,eta,mu,zet,s = subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);

    return xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp
end

function subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d)
    """
    This function subsolv solves the MMA subproblem:

    minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z +
             + SUM[ ci*yi + 0.5*di*(yi)^2 ],    
    subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
               alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.

    Input:  m, n, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d.
    Output: xmma,ymma,zmma, slack variables and Lagrange multiplers.
    """

    een = ones(n,1);
    eem = ones(m,1);
    epsi = 1;
    epsvecn = epsi*een;
    epsvecm = epsi*eem;
    x = 0.5*(alfa+beta);
    y = eem;
    z = 1;
    lam = eem;
    xsi = een./(x-alfa);
    xsi = max.(xsi,een);
    eta = een./(beta-x);
    eta = max.(eta,een);
    mu  = max.(eem,0.5*c);
    zet = 1;
    s = eem;
    itera = 0;
    while epsi > epsimin
      epsvecn = epsi*een;
      epsvecm = epsi*eem;
      ux1 = upp-x;
      xl1 = x-low;
      ux2 = ux1.*ux1;
      xl2 = xl1.*xl1;
      uxinv1 = een./ux1;
      xlinv1 = een./xl1;
      plam = p0 + P'*lam ;
      qlam = q0 + Q'*lam ;
      gvec = P*uxinv1 + Q*xlinv1;
      dpsidx = plam./ux2 - qlam./xl2 ;
      rex = dpsidx - xsi + eta;
      rey = c .+ d.*y - mu - lam;
      rez = a0 - zet .- a'*lam;
      relam = gvec .- a*z - y + s - b;
      rexsi = xsi.*(x-alfa) - epsvecn;
      reeta = eta.*(beta-x) - epsvecn;
      remu = mu.*y - epsvecm;
      rezet = zet*z - epsi;
      res = lam.*s - epsvecm;
      residu1 = [rex' rey' rez]';
      residu2 = [relam' rexsi' reeta' remu' rezet res']';
      residu = [residu1' residu2']';
      residunorm = sqrt.(residu'*residu)[1,1];
      residumax = maximum(abs.(residu));
      ittt = 0;
      while (residumax > 0.9*epsi) && (ittt < 200)
        ittt=ittt + 1;
        itera=itera + 1;
        ux1 = upp-x;
        xl1 = x-low;
        ux2 = ux1.*ux1;
        xl2 = xl1.*xl1;
        ux3 = ux1.*ux2;
        xl3 = xl1.*xl2;
        uxinv1 = een./ux1;
        xlinv1 = een./xl1;
        uxinv2 = een./ux2;
        xlinv2 = een./xl2;
        plam = p0 + P'*lam ;
        qlam = q0 + Q'*lam ;
        gvec = P*uxinv1 + Q*xlinv1;
        GG = P*sparse(Diagonal(uxinv2[:,1])) -Q*sparse(Diagonal(xlinv2[:,1]))#P*spdiags(uxinv2,0,n,n) - Q*spdiags(xlinv2,0,n,n);
        dpsidx = plam./ux2 - qlam./xl2 ;
        delx = dpsidx - epsvecn./(x-alfa) + epsvecn./(beta-x);
        dely = c .+ d.*y - lam - epsvecm./y;
        delz = a0 .- a'*lam .- epsi/z;
        dellam = gvec .- a*z - y - b .+ epsvecm./lam;
        diagx = plam./ux3 + qlam./xl3;
        diagx = 2*diagx + xsi./(x-alfa) + eta./(beta-x);
        diagxinv = een./diagx;
        diagy = d .+ mu./y;
        diagyinv = eem./diagy;
        diaglam = s./lam;
        diaglamyi = diaglam+diagyinv;
        if m < n
          blam = dellam .+ dely./diagy - GG*(delx./diagx);
          bb = [blam' delz]';
          Alam = sparse(Diagonal(diaglamyi[:,1]))+GG*sparse(Diagonal(diagxinv[:,1]))*GG'#spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
          AA = [Alam     a
                a'    -zet/z ];
          solut = AA\bb;
          dlam = solut[1:m];
          dz = solut[m+1];
          dx = -delx./diagx - (GG'*dlam)./diagx;
        else
          diaglamyiinv = eem./diaglamyi;
          dellamyi = dellam + dely./diagy;
          Axx = saprse(Diagonal(diagx[:,1])) + GG'*sparse(Diagonal(diaglamyiinv[:,1]))*GG#spdiags(diagx,0,n,n) + GG'*spdiags(diaglamyiinv,0,m,m)*GG;
          azz = zet/z + a'*(a./diaglamyi);
          axz = -GG'*(a./diaglamyi);
          bx = delx .+ GG'*(dellamyi./diaglamyi);
          bz  = delz .- a'*(dellamyi./diaglamyi);
          AA = [Axx   axz
                axz'  azz ];
          bb = [-bx' -bz]';
          solut = AA\bb;
          dx  = solut[1:n];
          dz = solut[n+1];
          dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
        end
    
        dy = -dely./diagy + dlam./diagy;
        dxsi = -xsi .+ epsvecn./(x-alfa) - (xsi.*dx)./(x-alfa);
        deta = -eta .+ epsvecn./(beta-x) + (eta.*dx)./(beta-x);
        dmu  = -mu .+ epsvecm./y - (mu.*dy)./y;
        dzet = -zet .+ epsi/z - zet*dz/z;
        ds   = -s .+ epsvecm./lam - (s.*dlam)./lam;
        xx  = [ y'  z  lam'  xsi'  eta'  mu'  zet  s']';
        dxx = [dy' dz dlam' dxsi' deta' dmu' dzet ds']';
        
        stepxx = -1.01*dxx./xx;
        stmxx  = maximum(stepxx);
        stepalfa = -1.01*dx./(x-alfa);
        stmalfa = maximum(stepalfa);
        stepbeta = 1.01*dx./(beta-x);
        stmbeta = maximum(stepbeta);
        stmalbe  = max.(stmalfa,stmbeta);
        stmalbexx = max.(stmalbe,stmxx);
        stminv = max.(stmalbexx,1);
        steg = 1/stminv;
    
        xold   =   x;
        yold   =   y;
        zold   =   z;
        lamold =  lam;
        xsiold =  xsi;
        etaold =  eta;
        muold  =  mu;
        zetold =  zet;
        sold   =   s;
    
        itto = 0;
        resinew = 2*residunorm;
        while resinew > residunorm && itto < 50
            itto = itto+1;
            x   =   xold + steg*dx;
            y   =   yold + steg*dy;
            z   =   zold + steg*dz;
            lam = lamold + steg*dlam;
            xsi = xsiold + steg*dxsi;
            eta = etaold + steg*deta;
            mu  = muold  + steg*dmu;
            zet = zetold + steg*dzet;
            s   =   sold + steg*ds;
            ux1 = upp-x;
            xl1 = x-low;
            ux2 = ux1.*ux1;
            xl2 = xl1.*xl1;
            uxinv1 = een./ux1;
            xlinv1 = een./xl1;
            plam = p0 + P'*lam ;
            qlam = q0 + Q'*lam ;
            gvec = P*uxinv1 + Q*xlinv1;
            dpsidx = plam./ux2 - qlam./xl2 ;
            rex = dpsidx - xsi + eta;
            rey = c .+ d.*y - mu - lam;
            rez = a0 - zet .- a'*lam;
            relam = gvec .- a*z - y + s - b;
            rexsi = xsi.*(x-alfa) - epsvecn;
            reeta = eta.*(beta-x) - epsvecn;
            remu = mu.*y - epsvecm;
            rezet = zet*z - epsi;
            res = lam.*s - epsvecm;
            residu1 = [rex' rey' rez]';
            residu2 = [relam' rexsi' reeta' remu' rezet res']';
            residu = [residu1' residu2']';
            resinew = sqrt.(residu'*residu)[1,1];
            steg = steg/2;
        end
      residunorm=resinew;
      residumax = maximum(abs.(residu));
      steg = 2*steg;
      end
      if ittt > 198
        epsi
        ittt
      end
    epsi = 0.1*epsi;
    end
    xmma   =   x;
    ymma   =   y;
    zmma   =   z;
    lamma =  lam;
    xsimma =  xsi;
    etamma =  eta;
    mumma  =  mu;
    zetmma =  zet;
    smma   =   s;
    return xmma,ymma,zmma,lamma,xsimma,etamma,mumma,zetmma,smma
end

function kktcheck(m,n,x,y,z,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval,dfdx,a0,a,c,d)
    """
    %  The left hand sides of the KKT conditions for the following
    %  nonlinear programming problem are calculated.
    %         
    %      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
    %    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
    %                xmax_j <= x_j <= xmin_j,    j = 1,...,n
    %                z >= 0,   y_i >= 0,         i = 1,...,m
    %*** INPUT:
    %
    %   m    = The number of general constraints.
    %   n    = The number of variables x_j.
    %   x    = Current values of the n variables x_j.
    %   y    = Current values of the m variables y_i.
    %   z    = Current value of the single variable z.
    %  lam   = Lagrange multipliers for the m general constraints.
    %  xsi   = Lagrange multipliers for the n constraints xmin_j - x_j <= 0.
    %  eta   = Lagrange multipliers for the n constraints x_j - xmax_j <= 0.
    %   mu   = Lagrange multipliers for the m constraints -y_i <= 0.
    %  zet   = Lagrange multiplier for the single constraint -z <= 0.
    %   s    = Slack variables for the m general constraints.
    %  xmin  = Lower bounds for the variables x_j.
    %  xmax  = Upper bounds for the variables x_j.
    %  df0dx = Vector with the derivatives of the objective function f_0
    %          with respect to the variables x_j, calculated at x.
    %  fval  = Vector with the values of the constraint functions f_i,
    %          calculated at x.
    %  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
    %          f_i with respect to the variables x_j, calculated at x.
    %          dfdx(i,j) = the derivative of f_i with respect to x_j.
    %   a0   = The constants a_0 in the term a_0*z.
    %   a    = Vector with the constants a_i in the terms a_i*z.
    %   c    = Vector with the constants c_i in the terms c_i*y_i.
    %   d    = Vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
    %     
    %*** OUTPUT:
    %
    % residu     = the residual vector for the KKT conditions.
    % residunorm = sqrt(residu'*residu).
    % residumax  = max(abs(residu)).
    """
        rex   = df0dx .+ dfdx .* lam - xsi + eta;
        rey   = c .+ d.*y - mu - lam;
        rez   = a0 - zet .- a'*lam;
        relam = fval .- a*z - y + s;
        rexsi = xsi.*(x-xmin);
        reeta = eta.*(xmax-x);
        remu  = mu.*y;
        rezet = zet*z;
        res   = lam.*s;
    
        residu1 = [rex' rey' rez]';
        residu2 = [relam' rexsi' reeta' remu' rezet res']';
        residu = [residu1' residu2']';
        residunorm = sqrt(residu'*residu);
        residumax = maximum(abs.(residu));
    return residu,residunorm,residumax
end
         