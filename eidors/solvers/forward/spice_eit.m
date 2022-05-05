function [Vn, In, nn] = spice_eit(netlist, freq)
% Vn = spice_eit(netlist, [freq])
% Modified Nodal Analysis (MNA) circuit solver
% for (almost/simple) spice netlists
% input
%   netlist - by cell array { [ part nodes...  args ], ... }
%             where parts are spice types:
%               "Rn" resistor; 2 nodes; args(1) = impedance
%               "Ln" inductor; 2 nodes; args(1) = inductance
%               "Cn" capacitor; 2 nodes; args(1) = capacitance
%               "Vn" voltage source; 2 nodes (-,+); args(1) = voltage
%               "In" current source; 2 nodes (-,+); args(1) = current
%               "Hn" current controlled voltage source, 4 nodes (-,+); args(1) = Vsrc, args(2) = gain
%               "En" voltage controlled voltage source, 4 nodes (-,+,d-,d+); args(1) = gain
%               "Fn" current controlled current source, 2 nodes (-,+); args(1) = Vsrc, args(2) = gain
%               "Gn" voltage controlled current source, 4 nodes (-,+,d-,d+); args(1) = gain
%   freq    - simulation frequency (default: 0 Hz, DC)
%   note: ground node is the lowest numbered node
% output
%   Vn      - nodal voltages by row [node voltage]
%   nn      - node numbers
%
% TODO reverse nodes for spice netlist (-,+) --> (+,-)
% TODO support DC voltage/current: [[DC] {value}]; currently supports {value}
% TODO support AC voltage/current: [AC {mag} [{phase}]]
%
% CITATION_REQUEST:
% AUTHOR: A Boyle and A Adler
% TITLE: Integrating Circuit Simulation with EIT FEM Models 
% JOURNAL: 19th International Conference on Biomedical Applications of Electrical Impedance Tomography, Edinburgh, UK
% YEAR: 2018
%

%  (C) 2017, 2018 A. Boyle, License: GPL version 2 or version 3

% based on http://matlabbyexamples.blogspot.ca/2011/11/circuit-solver-using-matlab-programming.html
% and https://www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA5.html
% and https://cseweb.ucsd.edu/classes/sp08/cse245/slides/let1.ppt
% and http://www.ecircuitcenter.com/SPICEsummary.htm
% and http://users.ecs.soton.ac.uk/mz/CctSim/chap1_4.htm
% and https://github.com/nik1106/MNA-MAT
%    (MNA-MAT is a different implementation of a spice-only linear simulator for DC/transient sims in matlab)
% and https://global.oup.com/us/companion.websites/fdscontent/uscompanion/us/static/companion.websites/9780199339136/Appendices/Appendix_B.pdf
% and http://www2.ensc.sfu.ca/~ljilja/papers/4C823d01.pdf
%    (non-linear components to ideal component circuits)

if ischar(netlist) && strcmp(netlist,'UNIT_TEST'); unittest(); return; end

if nargin < 2
   freq = 0;
end
w = 2*pi*freq;
s = j*w;

% collect nodes
assert(iscell(netlist), 'error: expected cell array for netlist');
nodes = zeros(length(netlist),2);
delM = zeros(length(netlist),1); % init: tracking L -> V converted inductors
M = 0; % number of voltage sources
P = 0;
for ii = 1:length(netlist)
   tt = netlist{ii};
   % save node numbers for later
   assert(isnumeric([tt{2:3}]), 'error: expected second and third field netlist to be node numbers');
   nodes(ii,:) = uint32([tt{2:3}]);
   id = tt{1};
   assert(ischar(id),'error: expected first field of netlist to be a component identifier');
   if abs(s) == 0
      % convert inductors to shorts at f=0; avoids 'matrix singular' errors
      switch id(1)
      case 'L'
         id = ['Vshort_f0_' id];
         tt{1} = id; tt{4} = 0;
         netlist{ii,1} = tt;
         delM(M+1) = 1;
      end
   end
   % count voltage sources and components;
   % this is a count of how many currents we will calculate
   if any(strcmp(id(1), {'V','E','H'}))
      M = M +1;
      Vidx{M} = id;
   end
   % extra entries in the sparse indices for large component stamps (dep. sources)
   switch id(1)
   case 'H'
      P = P+1;
   case 'E'
      P = P+2;
   end
end
nodes = uint32(nodes);
gnd = min(nodes(:));
delM = find(delM);

% remap node#s to be valid array indices
[ii] = unique(sort(nodes(:)));
jj = [1:length(ii)]';
mn = min(nodes(:));
map(ii-mn+1) = jj;
rmap(jj) = ii;

nodes = map(nodes - mn +1);
nodes = nodes -1; % gnd -> node#0


% save node numbers for output
nn = rmap(2:end)';


% build matrices
L = size(nodes,1);
N = length(nn); % number of nodes

%print_netlist('exec',netlist);

Aii = ones(L*4+P,1); Ajj=Aii; Avv=Aii*0; % init
Bii = ones(L*2,1); Bjj=Bii; Bvv=Bii*0; % init
Mn = N+1;
Pn = L*4+1;
for tt = 1:L
   row = netlist{tt}; val = row{4}; id = row{1};
   n1 = nodes(tt,1); n2 = nodes(tt,2);
   switch(id(1))
   case 'V' % independent voltage source
      idx = [ ((tt-1)*4+1):((tt-1)*4+4) ];
      % no off-diagonal entries for nodes connected to gnd
      vp = +1; vn = -1;
      if n1 == 0
         n1 = 1;
         vn = 0;
      elseif n2 == 0
         n2 = 1;
         vp = 0;
      end
      Aii(idx) = [ Mn Mn n1 n2 ];
      Ajj(idx) = [ n1 n2 Mn Mn ];
      Avv(idx) = [ vn vp vn vp ];
      idx = (tt-1)*2+1;
      Bii(idx) = [ Mn  ];
      Bvv(idx) = [ val ];
      Mn = Mn +1;
      continue;
   case 'F' % current controlled current source (CCCS)
      % from http://users.ecs.soton.ac.uk/mz/CctSim/chap1_4.htm
      Vsrc = row{4};
      Mi = find(strcmp(Vsrc, Vidx)) + N; % index of branch current
      assert(length(Mi) == 1,sprintf('failed to find %s for %s',Vsrc,id));
      gain = row{5}; assert(isnumeric(gain)); % gain
      idx = [ ((tt-1)*4+1):((tt-1)*4+2) ];
      % nodes connected to gnd?
      d1 = 1; d2 = 1;
      if n1 == 0
         n1 = 1;
         d1 = 0;
      end
      if n2 == 0
         n2 = 1;
         d2 = 0;
      end
      Aii(idx) = [  n1  n2 ];
      Ajj(idx) = [  Mi  Mi ];
      Avv(idx) = [ +d1 -d2 ]*gain;
      continue;
   case 'G' % voltage controlled current source (VCCS)
      % from https://cseweb.ucsd.edu/classes/sp08/cse245/slides/let1.ppt
      n1 = uint32(row{2}); % n-
      n2 = uint32(row{3}); % n+
      n3 = uint32(row{4}); % diff-
      n4 = uint32(row{5}); % diff+
      gain = row{6}; assert(isnumeric(gain)); % gain
      idx = [ ((tt-1)*4+1):((tt-1)*4+4) ];
      % no off-diagonal entries for nodes connected to gnd
      d1 = 1; d2 = 1; d3 = 1; d4 = 1;
      if n1 == 0
         n1 = 1;
         d1 = 0;
      end
      if n2 == 0
         n2 = 1;
         d2 = 0;
      end
      if n3 == 0
         n3 = 1;
         d3 = 0;
      end
      if n4 == 0
         n4 = 1;
         d4 = 0;
      end
      % A = [G B;C D] --> same as V; except C is changed
      %            C  C  B  B
      Aii(idx) = [  n1     n2     n1     n2 ];
      Ajj(idx) = [  n4     n3     n3     n4 ];
      Avv(idx) = [ -d1*d4 -d2*d3 +d1*d3 +d2*d4 ]*gain;
      continue;
   case 'E' % voltage controlled voltage source (VCVS)
      % from http://users.ecs.soton.ac.uk/mz/CctSim/chap1_4.htm
      n1 = uint32(row{2}); % n-
      n2 = uint32(row{3}); % n+
      n3 = uint32(row{4}); % diff-
      n4 = uint32(row{5}); % diff+
      gain = row{6}; assert(isnumeric(gain)); % gain
      idx = [ ((tt-1)*4+1):((tt-1)*4+4) Pn Pn+1 ];
      % no off-diagonal entries for nodes connected to gnd
      d1 = 1; d2 = 1; d3 = 1; d4 = 1;
      if n1 == 0
         n1 = 1;
         d1 = 0;
      end
      if n2 == 0
         n2 = 1;
         d2 = 0;
      end
      if n3 == 0
         n3 = 1;
         d3 = 0;
      end
      if n4 == 0
         n4 = 1;
         d4 = 0;
      end
      % A = [G B;C D] --> same as V; except C is changed
      %            C  C  B  B
      G = gain;
      Aii(idx) = [  Mn  Mn  n1  n2  Mn    Mn   ];
      Ajj(idx) = [  n1  n2  Mn  Mn  n3    n4   ];
      Avv(idx) = [ -d1 +d2 -d1 +d2 +d3*G -d4*G ];
      Pn = Pn +2;
      Mn = Mn +1;
      continue;
   case 'H' % current controlled voltage source (CCVS)
      % from http://users.ecs.soton.ac.uk/mz/CctSim/chap1_4.htm
      Vsrc = row{4};
      Mi = find(strcmp(Vsrc, Vidx)) + N; % index of branch current
      assert(length(Mi) == 1,sprintf('failed to find %s for %s',Vsrc,id));
      gain = row{5}; assert(isnumeric(gain)); % gain
      idx = [ ((tt-1)*4+1):((tt-1)*4+4) Pn ];
      % no off-diagonal entries for nodes connected to gnd
      d1 = 1; d2 = 1; d3 = 1; d4 = 1;
      if n1 == 0
         n1 = 1;
         d1 = 0;
      end
      if n2 == 0
         n2 = 1;
         d2 = 0;
      end
      Aii(idx) = [  Mn  Mn  n1  n2  Mn   ];
      Ajj(idx) = [  n1  n2  Mn  Mn  Mi   ];
      Avv(idx) = [ -d1 +d2 -d1 +d2  gain ];
      Pn = Pn +1;
      Mn = Mn +1;
      continue;
   case 'I' % independent current source
      idx = [ ((tt-1)*2+1):((tt-1)*2+2) ];
      ip = val; in = -val;
      if n1 == 0
         n1 = 1;
         in = 0;
      elseif n2 == 0
         n2 = 1;
         ip = 0;
      end
      Bii(idx) = [ n1 n2 ];
      Bvv(idx) = [ in ip ];
      continue; 
   case 'R' % resistor
      y = 1/val;
   case 'C' % capacitor
      y = s*val;
   case 'L' % inductor
      y = 1/(s*val);
   otherwise
      error(['unhandled netlist component type: ' id]);
   end
   if n1 == n2; error(sprintf('bad netlist @ line#%d: n1(%d)==n2(%d)',tt,n1,n2)); end

   % passive elements
   yn = -y; y1 = y; y2 = y;
   % no off-diagonal entries for nodes connected to gnd
   if n1 == 0
      n1 = 1;
      yn = 0;
      y1 = 0;
   elseif n2 == 0
      n2 = 1;
      yn = 0;
      y2 = 0;
   end
   idx = [ ((tt-1)*4+1):((tt-1)*4+4) ];
   Aii(idx) = [n1 n2 n1 n2];
   Ajj(idx) = [n2 n1 n1 n2];
   Avv(idx) = [yn yn y1 y2];
   % ... sparse-ified
   %A(n1,n2) = A(n1,n2)-y;
   %A(n2,n1) = A(n2,n1)-y;
   %A(n1,n1) = A(n1,n1)+y;
   %A(n2,n2) = A(n2,n2)+y;
end
% delete ground node
%idx = find(ii==gnd); ii(idx) = []; jj(idx) = []; vv(idx) = [];
%idx = find(jj==gnd); ii(idx) = []; jj(idx) = []; vv(idx) = [];
% build matrix
A = sparse(Aii,Ajj,Avv,N+M,N+M);
B = sparse(Bii,Bjj,Bvv,N+M,1);

%AA=full(A)
%BB=full(B)

% should be sparse and symmetric! check we use the right solver
X = A\B;
Vn = full(X(1:N));
In = full(X(N+1:N+M));
In(delM) = [];

function unittest()
   pass = 1;
   netlist_rlc = {{'R1', 0, 1, 1e3},
                  {'L1', 0, 1, 1},
                  {'C1', 0, 1, 5e-6},
                  {'RC2', 2, 1, 1},
                  {'R2', 3, 2, 1},
                  {'V1', 0, 3, 1}};
   tol = eps*10;
   for ii = 1:9
      switch(ii)
      case 1
         desc = 'voltage divider';
         netlist = {{'R1', 1, 2, 1},
                    {'R2', 0, 2, 1},
                    {'V0', 0, 1, 2}};
         Ev = [2,1]';
         Ei = [-1]';
         Enn = [1 2]';
         freq = 0;
      case 2
         desc = 'resistor network; gnd = n#4';
         % from http://matlabbyexamples.blogspot.ca/2011/11/circuit-solver-using-matlab-programming.html
         netlist = {{'R1', 1, 2, 2},
                    {'R2', 1, 3, 4},
                    {'R3', 2, 3, 5.2},
                    {'R4', 3, 4, 6},
                    {'R5', 2, 4, 3},
                    {'V1', 4, 1, 6},
                    {'V0', 0, 4, 0}};
         Ev = [6,3.6,3.6,0]';
         Ei = [-1.8,0]';
         Enn = [1 2 3 4]';
         freq = 0;
      case 3
         desc = 'reactive elements, f=0';
         netlist = netlist_rlc;
         Ev = [0.0,0.5,1.0]';
         Ei = [-0.5]';
         Enn = [1 2 3]';
         freq = 0;
      case 4
         desc = 'reactive elements, f=1';
         netlist = netlist_rlc;
         Ev = [0.90655 + 0.28793i;
               0.95328 + 0.14397i;
               1];
         Ei = -0.04672 + 0.14397i;
         Enn = [1 2 3]';
         freq = 1;
         tol = 6e-6;
      case 5
         desc = 'reactive elements, f=10';
         netlist = netlist_rlc;
         Ev = [0.99704 + 0.03105i;
               0.99852 + 0.01552i;
               1.0];
         Ei = -0.00148 + 0.01552i;
         Enn = [1 2 3]';
         freq = 10;
      case 6
         desc = 'voltage controlled current source VCCS G';
         netlist = {{'V1', 0, 1, 1},
                    {'R1', 0, 1, 1},
                    {'RL', 0, 2, 1},
                    {'G1', 0, 2, 0, 1, 5}};
         Ev = [1,-5]';
         Ei = [-1]';
         Enn = [1 2]';
         freq = 0;
      case 7
         desc = 'current controlled current source CCCS F';
         netlist = {{'V1', 0, 1, 1},
                    {'R1', 0, 1, 1},
                    {'RL', 0, 2, 1},
                    {'F1', 0, 2, 'V1', 5}};
         Ev = [1,-5]';
         Ei = [-1]';
         Enn = [1 2]';
         freq = 0;
      case 8
         desc = 'voltage controlled voltage source VCVS E';
         netlist = {{'V1', 0, 1, 1},
                    {'R1', 0, 1, 1},
                    {'RL', 0, 2, 1},
                    {'E1', 0, 2, 0, 1, 5}};
         Ev = [1,+5]';
         Ei = [-1,-5]';
         Enn = [1 2]';
         freq = 0;
      case 9
         desc = 'current controlled voltage source CCVS H';
         netlist = {{'V1', 0, 1, 1},
                    {'R1', 0, 1, 1},
                    {'RL', 0, 2, 1},
                    {'H1', 0, 2, 'V1', 5}};
         Ev = [1,+5]';
         Ei = [-1,-5]';
         Enn = [1 2]';
         freq = 0;
      otherwise
         error('oops');
      end
      desc = [sprintf('#%d ',ii) desc];

      disp('-----------------------------------------');
      print_netlist(desc,netlist);
      [v,i,nn]=spice_eit(netlist,freq);
      pass=cmp(pass,'voltages',v,Ev,tol);
      pass=cmp(pass,'currents',i,Ei,tol);
      pass=cmp(pass,'node#s',nn,Enn);
   end
   disp('-----------------------------------------');
   if pass
      disp('overall: PASS');
   else
      disp('overall: FAIL');
   end

function print_netlist(desc,n)
   fprintf('NETLIST: %s\n',desc);
   for ii =1:length(n)
      nn = n{ii}; id = nn{1};
      if any(strcmp(id(1),{'G','E'}))
         assert(length(nn) == 6,sprintf('bad netlist, row#%d has %d elements',ii,length(nn)));
         fprintf('  %-5s (%2d,%2d)<-(%2d,%2d)  %0.3f\n',nn{1},nn{2},nn{3},nn{4},nn{5},nn{6});
      elseif any(strcmp(id(1),{'F','H'}))
         assert(length(nn) == 5,sprintf('bad netlist, row#%d has %d elements',ii,length(nn)));
         fprintf('  %-5s (%2d,%2d)<-(%-5s)  %0.3f\n',nn{1},nn{2},nn{3},nn{4},nn{5});
      else
         assert(length(nn) == 4,sprintf('bad netlist, row#%d has %d elements',ii,length(nn)));
         fprintf('  %-5s (%2d,%2d)           %0.3f\n',nn{1},nn{2},nn{3},nn{4});
      end
   end

function pass=cmp(pass,str,X,Y,tol)
   if nargin < 5
      tol = 0;
   end
   if any(size(X) ~= size(Y))
      fprintf('TEST: %-30s fail; size mismatch (%d,%d)!=(%d,%d)\n',str,size(X),size(Y));
      pass=0;
      return
   end
   if tol == 0
      if any(X ~= Y)
         [ X Y ]
         fprintf('TEST: %-30s fail (%g)\n',str,abs(max(X(:)-Y(:))));
         pass=0;
         return;
      end
   else
      err=abs(X-Y);
      err(err<tol)=[];
      if length(err) > 0
         [ X Y ]
         fprintf('TEST: %-30s fail (%g, tol=%g)\n',str,norm(err),tol);
         pass=0;
         return;
      end
   end
   fprintf('TEST: %-30s pass\n',str);
