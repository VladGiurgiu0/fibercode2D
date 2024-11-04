function vorticity = compute_vorticity(gradTensor)
    % COMPUTE_VORTICITY computes the vorticity from the velocity gradient tensor
    [rows, cols] = size(gradTensor);
    vorticity = NaN(rows, cols);

    for i = 1:rows
        for j = 1:cols
            J = gradTensor{i, j};
            vorticity(i, j) = J(2, 1) - J(1, 2);
        end
    end
end