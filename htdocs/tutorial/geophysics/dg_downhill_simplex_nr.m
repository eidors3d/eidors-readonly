function [simplex y iter] =  dg_downhill_simplex_nr(func,simplex,y,ftol,itmax)
%
% simplex = [m1; m2; ...; mN]
%
    function amotry(fac)
        ndim= size(simplex,2);
        fac1= (1.0-fac)/ndim;
        fac2= fac1-fac;
        model_try= psum*fac1-simplex(ihi,:)*fac2;
        ytry= func(model_try);
        if ytry < y(ihi)
            y(ihi)= ytry;
            psum= psum-simplex(ihi,:)+model_try;
            simplex(ihi,:)= model_try(:);
        end
    end

ndim= size(simplex,2);
iter= 0;
psum= sum(simplex);
while 1==1
    [ylo,ilo]= min(y);
    [yhi,ihi]= max(y);
    ytmp= y(ihi);
    y(ihi)= y(ilo);
    [ynhi,inhi]= max(y);
    y(ihi)= ytmp;
    rtol= 2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+eps('double'));
    if rtol < ftol
        y([1 ilo])= y([ilo 1]);
        simplex([1 ilo])= simplex([ilo 1]);
        return;
    end
    if rem(iter,100) == 0
%         disp(['iter = ' num2str(iter)]);
    end
    if iter >= itmax
        disp('ITMAX exceeded in dg_downhill_simplex_nr');
        return;
    end
    amotry(-1.0);
    iter= iter+1;
    if ytry <= y(ilo)
        amotry(2.0);
        iter= iter+1;
    elseif ytry >= y(inhi)
        ysave= y(ihi);
        amotry(0.5);
        iter= iter+1;
        if ytry >= ysave
            simplex= 0.5*(simplex+repmat(simplex(ilo,:),size(simplex,1),1));
            for i= 1:ndim+1
                if i ~= ilo
                    y(i)= func(simplex(i,:));
                end
            end
            iter= iter+ndim;
            psum= sum(simplex);
        end
    end
end
end