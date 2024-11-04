function gradTensor = compute_velocity_gradient_tensor(V, U)
    % COMPUTE_VELOCITY_GRADIENT_TENSOR computes the velocity gradient tensor
    % for each point in the 2D velocity field.
    % U and V are the velocity components on a 2D grid.

    % Check if the sizes of U and V match
    if ~isequal(size(U), size(V))
        error('U and V must be the same size.');
    end

    % Compute the gradients of the velocity components
    [dUdy, dUdx] = gradient(U);
    [dVdy, dVdx] = gradient(V);

    % Initialize a cell array to store the velocity gradient tensor at each point
    gradTensor = cell(size(U));

    % Loop over each point in the grid to compute the velocity gradient tensor
    for i = 1:size(U, 1)
        for j = 1:size(U, 2)
            % Construct the velocity gradient tensor
            J = [dUdx(i,j), dVdx(i,j); dUdy(i,j), dVdy(i,j)];
            
            % Store the tensor in the cell array
            gradTensor{i,j} = J;
        end
    end
end
