function [xDB,yDB] = f_findOulet(surfgrid,Ac)
%==========================================================================

global parFAULT
%==========================================================================

Amask = surfgrid.A;
Amask.Z(Amask.Z<Ac) = NaN;

[ll,ii] = deg2utmsimple([parFAULT.Latl parFAULT.Latr],[parFAULT.Lonl parFAULT.Lonr],str2double(parFAULT.utm(1:2)),parFAULT.utm(3));

ptF_l = [find(ll(1)<surfgrid.xvec,1) find(ii(1)>surfgrid.yvec,1)];
ptF_r = [find(ll(2)<surfgrid.xvec,1) find(ii(2)>surfgrid.yvec,1)];
a     = (ptF_l(2)-ptF_r(2))/(ptF_l(1)-ptF_r(1));
b     = round(ptF_r(2) - a*ptF_r(1));

x = 2:surfgrid.A.size(2)-1;
y = ceil(b+a*x);

inc = 1;
for i = 1:length(x)
    
    m = max(max(Amask.Z(y(i)-1:y(i)+1,x(i)-1:x(i)+1)));
    
    if(isnan(m)==0)
        [ii,jj] = find(Amask.Z(y(i)-1:y(i)+1,x(i)-1:x(i)+1)==m);
        if(ii == 1); ii = -1; end
        if(ii == 2); ii = 0; end
        if(ii == 3); ii = 1; end
        if(jj == 1); jj = -1; end
        if(jj == 2); jj = 0; end
        if(jj == 3); jj = 1; end
        
        xDB(inc) = x(i)+jj;
        yDB(inc) = y(i)+ii;
        
        inc = inc + 1;
    end
end





