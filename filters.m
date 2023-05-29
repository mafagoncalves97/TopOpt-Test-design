% ======================================================================
%>@file filters.m
%>@brief Compute matrices and values for filtering operations
%>
%>@param Data (@b struct) Struct wiht the problem parameters
%>@param Data.nelx (@b integer) Number of elements in the horizontal direction
%>@param Data.nely (@b integer)  Number of elements in the vertical direction
%>@param Data.rmin (@b integer) Minimum radius for the element filter
%>@param Data.rnmin (@b integer) Minimum radius for the node filter
%>
%>@retval Filters (@b struct) Fields (H,Hs and Hn, Hns) Element and nodal
%matrices and values for filtering operations
%>
%>@details
% ======================================================================
function [Filters] = filters(Data)
rmin = Data.rmin; rnmin = Data.rnmin; nelx = Data.nelx; nely = Data.nely;

%%Prepare filter for elements
iH = ones(nelx* nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH); Hs = sum(H,2);

%%Prepare filter for nodals
inH = ones((nelx+1)*(nely+1)*(2*(ceil(rnmin)+1))^2,1);
jnH = ones(size(inH)); snH = zeros(size(inH)); k =0;
[elex,eley] = meshgrid(1.5:nelx+0.5,1.5:nely+0.5);
for in1 = 1:nelx+1
    for jn1 = 1:nely+1
        en1 = (in1-1)*(nely+1)+jn1;
        for in2 = max(in1-ceil(rnmin),1):min(in1+ceil(rnmin)-1,nelx)
            for jn2 = max(jn1-ceil(rnmin),1):min(jn1+ceil(rnmin)-1,nely)
                en2 = (in2-1)*nely+jn2; 
                k = k+1; 
                inH(k) = en1;
                jnH(k) = en2;
                snH(k) = max(0,rnmin-sqrt((in1-elex(jn2,in2))^2+(jn1-eley(jn2,in2))^2));
            end
        end
    end
end
Hn = sparse(inH,jnH,snH); Hns = sum(Hn,2);

%% Output filters
Filters.H = H;Filters.Hs = Hs;Filters.Hn = Hn;Filters.Hns = Hns;

end
%======================================================================
%>@file filters.m
%>@brief Compute matrices and values for filtering operations
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================