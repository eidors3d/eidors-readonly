ii = 1;
Y = calc_system_mat(img(ii)); Y= Y.E;
nn= num_nodes(img(ii));
A = Y(1:nn,1:nn); B= Y(1:nn,nn+1:end); D = Y(nn+1:end,nn+1:end);
Dprime  = D - B'*inv(A)*B;
RR = -1./(Dprime - tril(Dprime));
RR(isinf(RR)) = 0;
