function TVop=TV_operator_2D(msh);
% 03/12/00 By Andrea Borsic
%
% This function constructs Total Variation operator
%
% function TVop=tv_op(msh);

% $Id: TV_operator_2D.m,v 1.2 2008-03-16 00:57:50 aadler Exp $


num_tri=length(msh.TC);
num_nod=length(msh.PC);

list=[];

%We build a selection array, to index the T matrix

SEL=[ 1 2 3; 2 3 1];

INC=spalloc(num_nod,num_nod,num_tri*3);

for k=1:num_tri
   
   for j=1:3
      
      m=msh.TC(SEL(1,j),k);
      n=msh.TC(SEL(2,j),k);
      
      if m<n % ( only the upper triangular )
         
         if (INC(m,n)==0) INC(m,n)=k;
         else
            len=sqrt((msh.PC(1,m)-msh.PC(1,n))^2+(msh.PC(2,m)-msh.PC(2,n))^2);
				list(size(list,1)+1,:)=[len,k,INC(m,n)]; % Triangles and length of the shared edge are written into the list
         end % if then else
         
      else
         
         if (INC(n,m)==0) INC(n,m)=k;
         else
            len=sqrt((msh.PC(1,m)-msh.PC(1,n))^2+(msh.PC(2,m)-msh.PC(2,n))^2);
				list(size(list,1)+1,:)=[len,k,INC(n,m)]; % Triangles and length of the shared edge are written into the list
         end % if then else
         
         
      end % if m<n
      
   end % for j
   
end % for k

TVop=spalloc(num_tri,num_tri,2*length(list));

for i=1:length(list)
   
   TVop(i,list(i,2))=list(i,1);
   
   TVop(i,list(i,3))=-list(i,1);
    
end % for

TVl=list(:,1);

% Bye Bye

