function p = order_loop(pp,clk)
% takes a list of points and orders them clockwise, starting from the first

if nargin == 1
   clk = 1;
end

D = distmat(pp) + diag(inf*ones(length(pp),1));
p = zeros(size(pp));
p(1,:) = pp(1,:);
idx1 = 1;
for i = 1:length(pp)-1
   [jnk,idx2] = min(D(idx1,:));
   p(i+1,:) = pp(idx2,:);
   D(:,idx1) = [];
   D(idx1,:) = [];
   pp(idx1,:) = [];
   if idx2 > idx1
      idx1 = idx2-1;
   else
      idx1 = idx2;
   end
end
ctr = mean(p);
tmp = p - repmat(ctr,length(p),1);
th = unwrap(cart2pol(tmp(:,1),tmp(:,2)));

if th(1) < th(end)
   %counter-clockwise
   p = flipud(p);
end

if clk == -1
   p = flipud(p);
end
