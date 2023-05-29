% ======================================================================
%>@file lk.m
%>@brief Compute elemental stiffness matrix
%>
%>@param E0 (@b integer) Number of elements in the horizontal direction
%>@param v (@b integer)  Number of elements in the vertical direction
%>@param L (@b integer) Minimum radius for the element filter
%>@param A0 (@b integer) Minimum radius for the node filter
%>@param integration_type (@b integer) Minimum radius for the node filter
%>
%>@retval Ke (@b sparse matrix) Element and nodal matrix for filtering
%>@retval De (@b integer) Element and nodal values for filtering
%>@retval Be (@b integer) Element and nodal values for filtering
%>
%>@details
% ======================================================================
function [Ke, De, Be]=lk(Data)
% Returns the element stiffness matrix, shape funnction derivative matrix
% and elasticity tensor based on the integration type
 
v = Data.poisson; E0 = Data.E0;


%PLANE STRESS ELASTICITY MATRIX

D11 = (E0/(1-v^2));
D12 = (E0/(1-v^2)) * v;
D22 = (E0/(1-v^2));
D33 = (E0/(1-v^2))*(1-v)/2;
G12 = D33;

De = [D11 D12 0;
      D12 D22 0
      0    0  D33];
  

%SHAPE FUNCTION DERIVATIVE MATRIX

    if Data.integration == 'r' %reduced
        
        % AVERAGE OF SHAPE FUNCTION DERIVATIVE MATRIX OF INTEGRATION POINTS
        Be = (1/2/L)*[-1 0 1 0 1 0 -1 0; 0 -1 0 -1 0 1 0 1; -1 -1 -1 1 1 1 1 -1];
        Ke = Be'*De*Be*A0;
    
    elseif Data.integration == 'f' %full
        
        % SHAPE FUNCTION DERIVATIVE MATRIX OF 4 INTEGRATION POINTS 
        Be=zeros(3,8,4);
        a = 0.7895; b = 0.2105;
        Be(:,:,1) =[-a   0    a    0    b    0   -b    0;
                     0  -a    0   -b    0    b    0    a;
                    -a  -a   -b    a    b    b    a   -b];


        Be(:,:,2) =[-a   0    a    0    b    0   -b    0;
                     0  -b    0   -a    0    a    0    b;
                    -b  -a   -a    a    a    b    b   -b];


        Be(:,:,3) =[-b    0    b    0    a    0   -a    0;
                     0   -a    0   -b    0    b    0    a;
                    -a   -b   -b    b    b    a    a   -a];


        Be(:,:,4) = [-b    0    b    0    a    0   -a    0;
                      0   -b    0   -a    0    a    0    b;
                     -b   -b   -a    b    a    a    b   -a];

        Be=Be/Data.L;
  

        % PRECOMPUTED STIFFNESS MATRIX - PLANE STRESS - FULL INTEGRATION -
        % LINEAR
        Ke=zeros(8);
        Ke(1,2)= (D12+G12)/4;Ke(1,3)=-D11/3+G12/6;Ke(1,4)= (D12-G12)/4;
        Ke(1,5)=-(D11+G12)/6;Ke(1,6)=-(D12+G12)/4;Ke(1,7)= D11/6-G12/3;
        Ke(1,8)=-D12/4+G12/4;Ke(2,3)=-D12/4+G12/4;Ke(2,4)= D22/6-G12/3;
        Ke(2,5)=-(D12+G12)/4;Ke(2,6)=-(D22+G12)/6;Ke(2,7)= (D12-G12)/4;
        Ke(2,8)=-D22/3+G12/6;Ke(3,4)=-(D12+G12)/4;Ke(3,5)= D11/6-G12/3;
        Ke(3,6)= (D12-G12)/4;Ke(3,7)=-(D11+G12)/6;Ke(3,8)= (D12+G12)/4;
        Ke(4,5)=-D12/4+G12/4;Ke(4,6)=-D22/3+G12/6;Ke(4,7)= (D12+G12)/4;
        Ke(4,8)=-(D22+G12)/6;Ke(5,6)= (D12+G12)/4;Ke(5,7)=-D11/3+G12/6;
        Ke(5,8)= D12/4-G12/4;Ke(6,7)=-D12/4+G12/4;Ke(6,8)= D22/6-G12/3;
        Ke(7,8)=-D12/4-G12/4;
        Ke=Ke+Ke';
        Ke(1,1)=(D11+G12)/3;Ke(2,2)=(D22+G12)/3;Ke(3,3)=(D11+G12)/3;
        Ke(4,4)=(D22+G12)/3;Ke(5,5)=(D11+G12)/3;Ke(6,6)=D22/3+G12/3;
        Ke(7,7)=D11/3+G12/3;Ke(8,8)=D22/3+G12/3;

        Ke=Ke*(Data.L^2);
    end

end
%======================================================================
%>@file lk.m
%>@brief Compute element stiffness matrix
%>@details
%>
%>@author Mafalda Gonçalves
%>@date since 2020
%======================================================================