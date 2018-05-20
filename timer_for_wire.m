% Plots the time it takes to solve the General Hamiltonian up to an n size 
% qca wire. 

cell_num = 7;

sigmaZ = [1,0;0,-1];

P_d = linspace(-1, 1, 100);
y = 0.5;
times = zeros(1,cell_num);

for N = 1:cell_num
    tic
    polarization = zeros(100,1);
    for i = 1:100
        [V,D] = modularWire(N, y, P_d(i));
        polarization(i) = (V(:,1)')*paulizKron(sigmaZ,N,N)*V(:,1);
    end
    times(N) = toc;
end

plot(1:cell_num,times)
title('Calculation time as number of cells increases')
ylabel('Time')
xlabel('Number of cells in the wire.')


% plot(P_d,polarization)
% title(sprintf('Plot for input: N = %d,%d', N))
% xlabel('y')
% ylabel('Polarization')

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