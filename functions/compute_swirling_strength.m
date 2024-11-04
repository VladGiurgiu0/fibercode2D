function lambda_ci = compute_swirling_strength(gradTensor)
    % COMPUTE_SWIRLING_STRENGTH computes the swirling strength from the velocity gradient tensor
    [rows, cols] = size(gradTensor);
    lambda_ci = NaN(rows, cols);

    for i = 1:rows
        for j = 1:cols
            J = gradTensor{i, j};
            eigenvalues = eig(J);
            lambda_ci(i, j) = max(abs(imag(eigenvalues)));
        end
    end
end