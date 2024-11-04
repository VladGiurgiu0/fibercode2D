function [shear_rate] = compute_shear_rate(gradTensor)
    % COMPUTE_STRAIN_RATE computes the strain rate from the velocity gradient tensor
    [rows, cols] = size(gradTensor);
    shear_rate = NaN(rows, cols);

    for i = 1:rows
        for j = 1:cols
            J = gradTensor{i, j};

            %% V1
%             S = 0.5 * (J + J');
%             %strain_rate(i, j) = norm(S, 'fro');
% 
%             % Calculate the traceless part of A
%             D = S - (trace(S) / 2) * eye(2);
%             shear_rate(i, j) = norm(D, 'fro');

            %% V2
            shear_rate(i, j) = J(1,2) + J(2,1);
        end
    end
end
