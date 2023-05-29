% ======================================================================
%>@file materialBehavior.m
%>@brief Backward-Euler algorithm to compute stress and equivalent plastic
%>strain
%>@param estrain (@b integer) 
%>@param d_estrain (@b integer)  
%>@param D (@b integer) 
%>@param plastic_exist (@b integer) 
%>@param vxPhys (@b integer) 
%>@param penal (@b integer) 
%>@param Emin (@b integer) 
%>@param hv1 (@b integer) 
%>
%>@retval estress (@b array)
%>@retval Ct
%>@retval plastic_inc 
%>@retval ac (@b array)
%>@retval d_estrain_e (@b array) Element elastic strain increment
%>
%>@details
% ======================================================================
function [estress, Ct, plastic_inc, ac,d_estrain_e] = materialBehavior(estrain, d_estrain, D, plastic_exist,vxPhys, penal,Emin, hv1)
    %% Input parameters
    % estrain - actual element strain
    % d_estrain - element strain increment
    % E0, v - material properties
    % K, eps_0, sigma_0, n - Swift law parameters

    %Swift Law
    Ks = 979.46;  n = 0.194; sigma_0 = 355.0482; K_min = 300;
%     Ks = 529.5; n = 0.268; sigma_0 = 123.6;
    eps_0=(sigma_0/Ks)^(1/n);
    delta_lambda=0;
    if vxPhys>0.5        
        Ks = ((1-vxPhys)*Emin.^(penal-0.25)+vxPhys) * Ks;
    else
        Ks = 0.5 * Ks;
    end
    %% Start process

    %Elastic trial stress
    sigma_b = D * estrain + D * d_estrain;
    

    %Trial yield function (von Mises yield criteria)
    sigma_vm = (sigma_b(1)^2 + sigma_b(2)^2 - sigma_b(1)*sigma_b(2) + 3*sigma_b(3)^2)^(1/2);

    %Determine if actively yielding
    delta_lambda = plastic_exist;
    sigma_0_actual = Ks * (eps_0+ delta_lambda)^n;
    f = sigma_vm - sigma_0_actual;


    %Determine elastic and plastic strain and stress increments
    if f <= 0
        delta_lambda = 0;
        estress = sigma_b;
        d_estrain_e = d_estrain;
        Ct = D;
        plastic_inc=0;
        ac=zeros(3,1);

    elseif f > 0
        A = [2, -1, 0;-1, 2, 0;0, 0, 6];
        H = (Ks * n * (eps_0 + delta_lambda)^(n-1));
        a = [2 * sigma_b(1) - sigma_b(2);2 * sigma_b(2) - sigma_b(1); 6 * sigma_b(3)];
        a= a./(2*sigma_vm);
%         Al=(H*sigma_b'*a)/(sigma_0_actual);
        delta_lambda=f/(a'*D*a+H);
        sigma_c = sigma_b - delta_lambda .* (D * a);
%         sigma_0_actual = Ks * (eps_0+delta_lambda)^n;
        H = (Ks * n * (eps_0 + delta_lambda)^(n-1));
        sigma_vm = (sigma_c(1)^2 + sigma_c(2)^2 - sigma_c(1)*sigma_c(2) + 3*sigma_c(3)^2)^(1/2);
        f = sigma_vm - sigma_0_actual;
        
        i=1;

        %While the point is not in the yield locus
        while f > 1E-10
                        
            ac = [2 * sigma_c(1) - sigma_c(2);2 * sigma_c(2) - sigma_c(1);6 * sigma_c(3)];
            ac = ac./(2*sigma_vm);
            r = sigma_c - (sigma_b - delta_lambda * (D * ac));
            
            da_dsigma = (1/(2*sigma_vm))*A - (1/sigma_vm)*(ac*ac');
            Q = (eye(3) + delta_lambda*(D*da_dsigma));
           

            d_delta_lambda = (f - (((ac'*(Q\r)))))/((((ac'*(Q\D)))*ac) + H);
            d_sigma = -(Q\r + d_delta_lambda*(Q\D)*ac);
            delta_lambda = delta_lambda + d_delta_lambda;
            
%             sigma_0_actual = Ks * (eps_0+delta_lambda)^n;
            H = (Ks * n * (eps_0 + delta_lambda)^(n-1));
        
            
            sigma_c = sigma_c + d_sigma;
            sigma_vm = (sigma_c(1)^2 + sigma_c(2)^2 - sigma_c(1)*sigma_c(2) + 3*sigma_c(3)^2)^(1/2);
            f = sigma_vm - sigma_0_actual;

            
            i=i+1;
        end
       
        %Final calculations - elastic strain and stress
        ac = [2 * sigma_c(1) - sigma_c(2);2 * sigma_c(2) - sigma_c(1);6 * sigma_c(3)];
        ac = ac./(2*sigma_vm);
        d_estrain_e = d_estrain - delta_lambda*ac;

        estress = sigma_c;
        plastic_inc = delta_lambda;
        
        fprintf('Iterations: %d Stress: %.2f %.2f %.2f   VonMises: %.2f   Yield stress: %.2f  X: %.2f  Element: %d   Plastic: %.5f\n', i,sigma_c(1),sigma_c(2),sigma_c(3),sigma_vm, sigma_0_actual, vxPhys, hv1,plastic_inc)
                
        %Consistent tangent modular matrix
        R = Q\D;        
        H = H*(1-(plastic_exist+plastic_inc)*H);
        Ct = R*(eye(3) - ((ac*ac'*R')/(ac'*R*ac))+H);
    end

end
%======================================================================
%>@file materialBehavior.m
%>@brief Backward-Euler algorithm to compute stress and equivalent plastic
%>strain
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================