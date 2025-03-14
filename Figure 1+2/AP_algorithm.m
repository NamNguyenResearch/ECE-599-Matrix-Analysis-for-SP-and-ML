function [a,a_Kronecker,iter]=AP_algorithm(h,iter_max,error_bound,Us,K)
   %====================== Time line of coding ===========================
    % Programmer: Nam Nguyen
    % Oregon State University, Corvallis, Oregon, United State
    % Date: May 25, 2022

    %====================== Usage =========================================

    % [a,a_Kronecker,iter]=AP_algorithm(h,iter_max,error_bound,Us,K)

    %====================== Parameters ====================================
    % Outputs: 
    % a             :   a column vector of the estimated mixing matrix
    % a_Kronecker   :   Kronecker product of column vector
    % iter          :   number of iterations
    %----------------------------------------------------------------------
    % Inputs:
    % h             :   a vector belongs to the range space of the basis matrix Us
    % iter_max      :   maximal number of iterations
    % error_bound   :   tolerance between the currently estimated channel and the previously estimeated channel
    % Us            :   Khatri-Rao subspace
    % K             :   number of sources

    
    %% Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    value_previous=inf;

    for i=1:iter_max
        H=reshape(h,K,K);               % inverse vectorization 
        H=(H+H')/2; 
        
        % Calculate eigen value and eigen vector of H %%%%%%%%%%%%%%%%%%%%%
        [eigvec,eigval]=eig(H);
        [eigval, order_eig]=sort(diag(abs(eigval)),'descend');
        eigvec=eigvec(:,order_eig);
        
        % Close-form solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a=eigvec(:,1);                  % unit-2-norm eigenvector
        a_Kronecker=kron(conj(a),a);    % Kronecker product
        alpha=eigval(1)./abs(eigval(1));

        h=alpha*(Us*Us')*a_Kronecker;   % solution of h
        
        % Current objective  value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        value_current=sum(abs(h-alpha*a_Kronecker).^2); 
        
        % Stopping criterion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i>1
            if abs(value_previous-value_current)/abs(value_previous)<error_bound
               iter= i; 
               return; 
            end
        end
        
        value_previous=value_current;
    end
    
    iter= iter_max;
end