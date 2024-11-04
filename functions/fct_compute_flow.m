function [] = fct_compute_flow(p,iiii)
%% load flow data
switch p.which_flow_field
    case 'pairs'
        parfor i=1:p.nr_frames
            data = load([p.main_folder,'/only_tracers/out_PaIRS/','out_',sprintf('%04.f',i-1),'.mat']);
            U_temp(:,:,i) = data.U;
            V_temp(:,:,i) = data.V;
            X_temp(:,:,i) = data.X;
            Y_temp(:,:,i) = data.Y;
        end

    case 'raft'
        file_name = [p.main_folder,'only_tracers/out_RAFT/','Dataset_RealPIV_TU_Wien_Re1000_v5_pre_proc_part_',num2str(iiii-1,'%03.f'),'_full.hdf5'];
        parfor i=1:p.nr_frames
            data = h5read(file_name,['/result_',num2str(i-1,'%05.f')]);
            U_IMG_temp(:,:,i) = data(:,:,1)';
            V_IMG_temp(:,:,i) = data(:,:,2)';
        end
end

%% filter flow_data in time
if p.filter_in_time==1
    switch p.which_flow_field
        case 'pairs'
                U_temp = smoothdata(U_temp,3,p.filter_time_type,p.filter_time_kernel);
                V_temp = smoothdata(V_temp,3,p.filter_time_type,p.filter_time_kernel);
        case 'raft'
                U_IMG_temp = smoothdata(U_IMG_temp,3,p.filter_time_type,p.filter_time_kernel);
                V_IMG_temp = smoothdata(V_IMG_temp,3,p.filter_time_type,p.filter_time_kernel);
    end
end


%% process flow field
switch p.which_flow_field
    case 'pairs'
        % image reference frame
        temp_img = imread([p.main_folder,'/only_tracers/','tracers_0001.tif']);
        x_img = 1:size(temp_img,2);
        y_img = 1:size(temp_img,1);

        [X_IMG, Y_IMG] = meshgrid(x_img,y_img);
        U_IMG_temp = []; V_IMG_temp = [];

        if p.interpolate_velocity==1
            parfor i=1:p.nr_frames
                U_IMG_temp(:,:,i) = interp2(X_temp(:,:,i),Y_temp(:,:,i),U_temp(:,:,i),X_IMG,Y_IMG,p.interp_type);
                V_IMG_temp(:,:,i) = interp2(X_temp(:,:,i),Y_temp(:,:,i),V_temp(:,:,i),X_IMG,Y_IMG,p.interp_type); 
            end
        elseif p.interpolate_velocity==0
            parfor i=1:p.nr_frames
                U_IMG_temp(:,:,i) = padarray(repelem(U_temp(:,:,i),p.IW_size,p.IW_size), [p.IW_size/2, p.IW_size/2], NaN, 'both');
                V_IMG_temp(:,:,i) = padarray(repelem(V_temp(:,:,i),p.IW_size,p.IW_size), [p.IW_size/2, p.IW_size/2], NaN, 'both');
            end
        end
    case 'raft'
end

% check
if p.plot_images==1
    figure(1); clf; set(gcf,'Position',[100 100 800 800])
    for i=1:p.nr_frames
        subplot(1,2,1)
        imagesc(U_temp(:,:,i)); daspect([1 1 1])
        
        subplot(1,2,2)
        imagesc(U_IMG_temp(:,:,i)); daspect([1 1 1])
        pause
    end
end

%% filter flow field
if p.filter_velocity_field ==1
    switch p.which_flow_field
        case 'pairs'
                if p.interpolate_velocity==1
                    parfor i=1:p.nr_frames
                        U_IMG_temp(:,:,i) = imgaussfilt(U_IMG_temp(:,:,i),p.filter_sigma);
                        V_IMG_temp(:,:,i) = imgaussfilt(V_IMG_temp(:,:,i),p.filter_sigma);
                    end
                elseif p.interpolate_velocity==0
                    parfor i=1:p.nr_frames
                        U_temp(:,:,i) = imgaussfilt(U_temp(:,:,i),p.filter_sigma);
                        V_temp(:,:,i) = imgaussfilt(V_temp(:,:,i),p.filter_sigma);
                        U_IMG_temp = [];
                        V_IMG_temp = [];
                        U_IMG_temp(:,:,i) = padarray(repelem(U_temp(:,:,i),p.IW_size,p.IW_size), [p.IW_size/2, p.IW_size/2], NaN, 'both');
                        V_IMG_temp(:,:,i) = padarray(repelem(V_temp(:,:,i),p.IW_size,p.IW_size), [p.IW_size/2, p.IW_size/2], NaN, 'both');
                    end
                end
        case 'raft'
            parfor i=1:p.nr_frames
                U_IMG_temp(:,:,i) = imgaussfilt(U_IMG_temp(:,:,i),p.filter_sigma);
                V_IMG_temp(:,:,i) = imgaussfilt(V_IMG_temp(:,:,i),p.filter_sigma);
            end
    end
end

%% compute vorticity, lambda_ci, strain rate
switch p.which_flow_field
    case 'pairs'
        if p.interpolate_velocity==0
            parfor i=1:p.nr_frames
                gradTensor = compute_velocity_gradient_tensor(U_temp(:,:,i), V_temp(:,:,i));
                vorticity_IMG_temp(:,:,i) = single(padarray(repelem(compute_vorticity(gradTensor),p.IW_size,p.IW_size), [p.IW_size/2, p.IW_size/2], NaN, 'both'));
                lambda_ci_IMG_temp(:,:,i) = single(padarray(repelem(compute_swirling_strength(gradTensor),p.IW_size,p.IW_size), [p.IW_size/2, p.IW_size/2], NaN, 'both'));
                %strain_rate_IMG_temp(:,:,i) = single(padarray(repelem(compute_strain_rate(gradTensor),p.IW_size,p.IW_size), [p.IW_size/2, p.IW_size/2], NaN, 'both'));
                shear_rate_IMG_temp(:,:,i) = single(padarray(repelem(compute_shear_rate(gradTensor),p.IW_size,p.IW_size), [p.IW_size/2, p.IW_size/2], NaN, 'both'));
                
                %disp(['Computed vorticity, lambda_ci, and strain_rate for frame number ',num2str(i), ' out of ',num2str(p.nr_frames)])
            end

        elseif p.interpolate_velocity==1
            parfor i=1:p.nr_frames
                gradTensor = compute_velocity_gradient_tensor(U_IMG_temp(:,:,i), V_IMG_temp(:,:,i));
                vorticity_IMG_temp(:,:,i) = single(compute_vorticity(gradTensor));
                lambda_ci_IMG_temp(:,:,i) = single(compute_swirling_strength(gradTensor));
                %strain_rate_IMG_temp(:,:,i) = single(compute_strain_rate(gradTensor));
                shear_rate_IMG_temp(:,:,i) = single(compute_shear_rate(gradTensor));

                %disp(['Computed vorticity, lambda_ci, and strain_rate for frame number ',num2str(i), ' out of ',num2str(p.nr_frames)])
            end

        end
    case 'raft'
        parfor i=1:p.nr_frames
            gradTensor = compute_velocity_gradient_tensor(U_IMG_temp(:,:,i), V_IMG_temp(:,:,i));
            vorticity_IMG_temp(:,:,i) = compute_vorticity(gradTensor);
            lambda_ci_IMG_temp(:,:,i) = compute_swirling_strength(gradTensor);
            %strain_rate_IMG_temp(:,:,i)= compute_strain_rate(gradTensor);
            shear_rate_IMG_temp(:,:,i)= compute_shear_rate(gradTensor);

            %disp(['Computed vorticity, lambda_ci, and strain_rate for frame number ',num2str(i), ' out of ',num2str(p.nr_frames)])
        end
end

%% save flow fields
for i=1:p.nr_frames
    switch p.which_flow_field
        case 'pairs'
            U_IMG = single(U_IMG_temp(:,:,i));
            V_IMG = single(V_IMG_temp(:,:,i));
            vorticity_IMG = vorticity_IMG_temp(:,:,i);
            lambda_ci_IMG = lambda_ci_IMG_temp(:,:,i);
            %strain_rate_IMG = strain_rate_IMG_temp(:,:,i);
            shear_rate_IMG = shear_rate_IMG_temp(:,:,i);

            save([p.main_folder,'only_tracers/out_PaIRS/','flow_',num2str(i-1,'%04d'),'.mat'], ...
                'U_IMG','V_IMG','vorticity_IMG','lambda_ci_IMG','shear_rate_IMG','p');

        case 'raft'
            U_IMG = single(U_IMG_temp(:,:,i));
            V_IMG = single(V_IMG_temp(:,:,i));
            vorticity_IMG = vorticity_IMG_temp(:,:,i);
            lambda_ci_IMG = lambda_ci_IMG_temp(:,:,i);
            %strain_rate_IMG = strain_rate_IMG_temp(:,:,i);
            shear_rate_IMG = shear_rate_IMG_temp(:,:,i);

            save([p.main_folder,'only_tracers/out_RAFT/','flow_',num2str(i-1,'%04d'),'.mat'], ...
                'U_IMG','V_IMG','vorticity_IMG','lambda_ci_IMG','shear_rate_IMG','p');
    end
end
end