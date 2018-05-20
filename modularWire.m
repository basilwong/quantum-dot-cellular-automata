%% The below function simulates a 'wire' of cells. 
% The length of the wire can be modified with the N input.
function [V,D] = modularWire(N,y,P_d)
    
    sigmaZ = [1,0;0,-1];
    sigmaX = [0,1;1,0];
    
    % Making the first term of the general hamiltonian.
    H_y = 0;
    for index1 = 1:N
        H_y_new = -1*y*paulizKron(sigmaX,index1,N);
        H_y = H_y + H_y_new;
    end
    
   % Making the second term of the general hamiltonian.
    H_P_d = 0;
    for index2 = 1:N
        H_P_d_new = -0.5*P_d*get_E_k(index2,0)*paulizKron(sigmaZ,index2,N);
        H_P_d = H_P_d + H_P_d_new;
    end
    
    % Making the third term of the general hamiltonian. 
    H_E_k = 0;
    for indexi = 1:N-1
        for indexj = indexi+1:N
            H_E_k_new = -0.5*get_E_k(indexj-indexi,0)*paulizKron(sigmaZ,indexi,N)*paulizKron(sigmaZ,indexj,N);
            H_E_k = H_E_k + H_E_k_new;
        end
    end
    
    H = H_y + H_P_d + H_E_k;
    
    [V,D] = eig(H);
    D = diag(D);
    
end 

%% This function finds E_k between cells. Takes cell distance as input.
% The approximate equation for E_k is cos(4*theta)/r^5.
function E_k = get_E_k(inline_distance, perpendiculat_distance)

    % Finding total distance between cells. 
    r = sqrt(inline_distance^2 + perpendiculat_distance^2);
    
    % Finding angle between cells. 
    theta = atan(perpendiculat_distance/inline_distance);
    
    % Plugging values into final equation.
    E_k = cos(4*theta)/r^5;
end

% Function for putting a matrix within a kroeneker multiplication of
% identitiy matrices. 
function matrixR = paulizKron(A,i,N)

    I = eye(2);
    
    if(i>2)
        matrixR = eye(2^(i-1));
        matrixR = kron(matrixR,A);
    elseif(i==1)
        matrixR = A;
    elseif(i==2)
        matrixR = kron(I,A);
        
    end
    
    for k=i:N
        if (k<N)
            matrixR = kron(matrixR,I);
        end
    end
end

%     % Initiating the Hamiltonian of the first cell.
%     E_k = 0;
%     for p = 1:N
%         E_k =   E_k + get_E_k(p,0);
%     end
%     E_k =   E_k + get_E_k(1,0);
%             
%     H = [-1*E_k/2, -1*y;
%         -1*y, E_k/2];
%     
%     E_k = 0;
%     for i = 2:N
%         if i ~= N
%             for j = 1:N-i
%                 E_k =   E_k + get_E_k(j,0);
%             end
%         end
%         if i ~= 1
%             for k = 1:i-1
%                 E_k =   E_k + get_E_k(k,0);
%             end
%         end
%         
%         % Add the influence of the driver.
%         E_k = E_k + P_d*get_E_k(i,0);
%         
%         H_new = [-1*E_k/2, -1*y;
%         -1*y, E_k/2];
%         
%         H = kron(H,H_new);
%     
%     end
%     [V,D] = eig(H);
%     D = diag(D);
% end
