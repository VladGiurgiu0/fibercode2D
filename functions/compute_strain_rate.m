function [strain_rate] = compute_strain_rate(gradTensor)
    % COMPUTE_STRAIN_RATE computes the strain rate from the velocity gradient tensor
    [rows, cols] = size(gradTensor);
    strain_rate = NaN(rows, cols);
    for i = 1:rows
        for j = 1:cols
            J = gradTensor{i, j};
            S = 0.5 * (J + J');
            strain_rate(i, j) = norm(S, 'fro');
        end
    end
end
