%% The below function generates the general Hamiltonian for a qca circuit 
% majority gate. 
function [V,D] = generalHamiltonian(y, P_D_1, P_D_2, P_D_3)    

    sigmaZ = [1,0;0,-1];
    sigmaX = [0,1;1,0];
    
% Making the first term of the general hamiltonian.
    H_y = 0;
    for index1 = 1:11
        H_y_new = -1*y*paulizKron(sigmaX,index1,11);
        H_y = H_y + H_y_new;
    end
    
% Making the second term of the general hamiltonian.
    H_P_d = 0;
    for index2 = 1:11
        H_P_d_new = -0.5*drivers_E_k(index2, P_D_1, P_D_2, P_D_3)*paulizKron(sigmaZ,index2,11);
        H_P_d = H_P_d + H_P_d_new;
    end
    
% Making the third term of the general hamiltonian. 
    H_E_k = 0;
    for indexi = 1:10
        for indexj = indexi+1:11
            H_E_k_new = -0.5*getKinkEnergy(indexi,indexj)*paulizKron(sigmaZ,indexi,11)*paulizKron(sigmaZ,indexj,11);
            H_E_k = H_E_k + H_E_k_new;
        end
    end
    
    H = H_y + H_P_d + H_E_k;
    
    [V,D] = eig(H);
    D = diag(D);
    
end

%% This function finds the kink energy for the interaction between cells. 
function kinkEnergy = getKinkEnergy(indexi,indexj)
    if indexi <= 6
        if indexj <= 6
            kinkEnergy = get_E_k(indexj - indexi,0);
        else
            if indexi <= 4
                if indexj == 11
                    kinkEnergy = get_E_k(5-indexi,1);
                else
                    kinkEnergy = get_E_k(5-indexi,11-indexj);
                end
            elseif indexi == 5
                if indexj == 11
                    kinkEnergy = get_E_k(1,0);
                else
                    kinkEnergy = get_E_k(11-indexj,0);
                end
            elseif indexi == 6
                if indexj == 11
                    kinkEnergy = get_E_k(1,1);
                else
                    kinkEnergy = get_E_k(11-indexj,1);
                end
            end
        end
    else
        kinkEnergy = get_E_k(indexj - indexi,0);
    end
end
%% This function determines the driver influence for the specfied majority 
% gate.
function sum_E_k = drivers_E_k(index,P_D_1, P_D_2, P_D_3)
     switch index
        case 1
            P_1_E_k = get_E_k(1,0)*P_D_1;
            P_2_E_k = get_E_k(4,5)*P_D_2;
            P_3_E_k = get_E_k(6,0)*P_D_3;
            
        case 2
            P_1_E_k = get_E_k(2,0)*P_D_1;
            P_2_E_k = get_E_k(3,5)*P_D_2;
            P_3_E_k = get_E_k(5,0)*P_D_3;
        
        case 3
            P_1_E_k = get_E_k(3,0)*P_D_1;
            P_2_E_k = get_E_k(2,5)*P_D_2;
            P_3_E_k = get_E_k(4,0)*P_D_3;
            
        case 4
            P_1_E_k = get_E_k(4,0)*P_D_1;
            P_2_E_k = get_E_k(1,5)*P_D_2;
            P_3_E_k = get_E_k(3,0)*P_D_3;
            
        case 5
            P_1_E_k = get_E_k(5,0)*P_D_1;
            P_2_E_k = get_E_k(0,5)*P_D_2;
            P_3_E_k = get_E_k(2,0)*P_D_3;
            
        case 6
            P_1_E_k = get_E_k(6,0)*P_D_1;
            P_2_E_k = get_E_k(1,5)*P_D_2;
            P_3_E_k = get_E_k(1,0)*P_D_3;
            
        case 7
            P_1_E_k = get_E_k(4,5)*P_D_1;
            P_2_E_k = get_E_k(1,0)*P_D_2;
            P_3_E_k = get_E_k(4,2)*P_D_3;
             
        case 8
            P_1_E_k = get_E_k(3,5)*P_D_1;
            P_2_E_k = get_E_k(2,0)*P_D_2;
            P_3_E_k = get_E_k(3,2)*P_D_3;
            
        case 9 
            P_1_E_k = get_E_k(2,5)*P_D_1;
            P_2_E_k = get_E_k(3,0)*P_D_2;
            P_3_E_k = get_E_k(2,2)*P_D_3;           
            
        case 10
            P_1_E_k = get_E_k(1,5)*P_D_1;
            P_2_E_k = get_E_k(4,0)*P_D_2;
            P_3_E_k = get_E_k(1,2)*P_D_3;
            
        case 11
            P_1_E_k = get_E_k(1,5)*P_D_1;
            P_2_E_k = get_E_k(6,0)*P_D_2;
            P_3_E_k = get_E_k(1,2)*P_D_3;
    end
          
    sum_E_k = P_1_E_k + P_2_E_k + P_3_E_k;
end

%% This function finds E_k between cells. Takes cell distance as input.
% The approximate equation for E_k is cos(4*theta)/r^5.
function E_k = get_E_k(inline_distance, perpendiculat_distance)

    % Finding total distance between cells. 
    r = sqrt(inline_distance^2 + perpendiculat_distance^2);
    
    % Finding angle between cells. 
    if inline_distance >= perpendiculat_distance
        theta = atan(perpendiculat_distance/inline_distance);
    else
        theta = atan(inline_distance/perpendiculat_distance);
    end
    
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
        
        
