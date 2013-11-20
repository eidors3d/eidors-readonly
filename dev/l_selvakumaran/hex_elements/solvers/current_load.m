function [I]=current_load(v)
%no_el- no. of electrodes
%patt- current injection pattern (0-opposite/1-adjacent)
%int- intensity of the current applied
no_el=v.n_elec;
 m_ind=[no_el,2];
 if v.inj== '{ad}'
     for i=1:no_el;
         m_ind(i,1)=i;
         if i ~=no_el;
             m_ind(i,2)=i+1;
         else
             m_ind(i,2)=1;
         end
     end
 end
 if v.inj== '{op}'
     for i=1:no_el;
         m_ind(i,1)=i;
     end
     for i=1:no_el/4;
         m_ind(i,2)=(3/4)*no_el+1-i;
     end
     for i=(no_el/4)+1: no_el/2;
         m_ind(i,2)=(no_el+1-(i-(no_el/4)));
     end
     for i=(no_el/2)+1:(3*no_el/4);
         m_ind(i,2)=(no_el/4)-i+(no_el/2)+1;   
     end
     for i=(3*no_el/4)+1:no_el;
         m_ind(i,2)=(no_el/2)-i+(3*no_el/4)+1;   
     end
 end
 %m_ind
 I=[];
 for i=1:no_el
     m_n=zeros(no_el,1);
     m_n(m_ind(i,1))=-1;
     m_n(m_ind(i,2))=1;
     I=[I,m_n];
 end
     

end
