% Using the function for generating the Hamiltonian, we check to see if the
% output is a true majority gate. 

sigmaZ = [1,0;0,-1];

% A = 1;
B = 1;
C = -1;

A = linspace(-1, 1, 100);
% B = linspace(-1, 1, 100);
% C = linspace(-1, 1, 100);
polarization = zeros(100,1);
for i = 1:100
    [V,D] = generalHamiltonian(0.01, A(i), B, C);
    polarization(i) = (V(:,1)')*paulizKron(sigmaZ,11,11)*V(:,1);
end

plot(A,polarization)
title(sprintf('Plot for input: %d,%d,%d', A(100),B,C))
xlabel('y')
ylabel('Polarization')

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