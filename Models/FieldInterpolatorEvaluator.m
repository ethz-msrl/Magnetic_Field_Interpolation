classdef FieldInterpolatorEvaluator < handle
    %FIELDINTERPOLATOREVALUATOR Used For Evaluating the Performance of a
    % magnetic field interpolation model
    %   Copyright 2020, ETH Zurich, Multi Scale Robotics Lab, Samuel Charreyron
    
    properties (SetAccess = protected)
        NodePositions
        NodeValues
        Model
        FieldsPredicted
    end
    
    methods
        
        function [obj] = FieldInterpolatorEvaluator(model, nodes, values)
            % Constructor
            %   Args:
            %       model (FieldInterpolator): the model to evaluate
            %       nodes (4D array): the positions at which to interpolate
            %       the dimensions should be X,Y,Z,3
            %       values (4D array): the field values to interpolate
            %       the dimensions should be X,Y,Z,3
            obj.NodePositions = nodes;
            obj.NodeValues = values;
            
            if ndims(nodes) ~= 4
                error('nodes should be a 4D array')
            end
            
            if ndims(values) ~= 4
                error('values should be a 4D array')
            end
            
            obj.Model = model;
            
            if numel(nodes) ~= numel(values)
                error('Node positions and values dont have the same number of elements');
            end

            obj.FieldsPredicted = model.getFieldsAtPositions(obj.NodePositions);
            
        end
        
        function setModel(obj, model)
            obj.Model = model;
            obj.FieldsPredicted = model.getFieldsAtPositions(obj.NodePositions);
        end
        
        function [r2] = get_r2(obj)
            % Gets the R2 score over the 3 dimensions
            % Returns:
            %   r2: a 3D array with an r2 score for each dimension
          
            nv = reshape(obj.NodeValues, [], 3);
            mu_nv = mean(nv, 1);
            ss_tot = sum((nv - mu_nv).^2, 1);
            ss_res = sum((nv - reshape(obj.FieldsPredicted, [], 3)).^2, 1);
            r2 = 1 - ss_res ./ ss_tot;
        end
        
        function [rmse] = get_rmse(obj)
            % Gets the RMSE score over the 3 dimensions
            % Returns:
            %   rmse: a 3D array with an rmse score in each dimension
            nv = reshape(obj.NodeValues, [], 3);
            nv_ = reshape(obj.FieldsPredicted, [], 3);
            rmse = sqrt(mean((nv - nv_).^2,1));
        end
        
        function [nmse] = get_nrmse(obj)
            % Gets the RMSE normalized by the maximum field magnitude
            
            nv = reshape(obj.NodeValues, [], 3);
            nv_ = reshape(obj.FieldsPredicted, [], 3);
            rmse = sqrt(mean((nv - nv_).^2,1));
            nmse = rmse ./ (max(nv_) - min(nv_));
        end
        
        function [mae] = get_mae(obj)
            % Gets the MAE score
            nv = reshape(obj.NodeValues, [], 3);
            nv_ = reshape(obj.FieldsPredicted, [], 3);
            mae = mean(abs(nv - nv_));
        end
        
        function [nmae] = get_nmae(obj)
            % Gets the MAE normalized by the maximum field magnitude
            nv = reshape(obj.NodeValues, [], 3);
            nv_ = reshape(obj.FieldsPredicted, [], 3);
            mae = mean(abs(nv - nv_));
            nmae = mae ./ (max(nv_) - min(nv_));
        end
        
        function [fh] = plot_quiver_xy(obj, varargin)
            % Plots a 2D quiver plot on a plane slice parallel to the xy
            % plane
            
            % Args:
            %   idx: (integer) the index of the xy plane in NodePositions
            %        must be between 1 and size(NodePositions,1)
            %   scale: (float) the scale of the quiver arrows
            % Returns:
            %   fh: the figure handle
            Args = {floor(size(obj.NodePositions,1)/2), 0.1};
            Args(1:(nargin-1)) = varargin;
            idx = Args{1};
            fh = figure;
            q1 = quiver(obj.NodePositions(:,:,idx,1), obj.NodePositions(:,:,idx,2),...
                obj.FieldsPredicted(:,:,idx,1), obj.FieldsPredicted(:,:,idx,2), 0);
            hold on;
            q2 = quiver(obj.NodePositions(:,:,idx,1), obj.NodePositions(:,:,idx,2),...
                obj.NodeValues(:,:,idx,1), obj.NodeValues(:,:,idx,2), 0);
            hold off;
            
            legend('Predicted', 'Reference');

            scale = Args{2};
            qU1 = get(q1, 'UData');
            qV1 = get(q1, 'VData');
            set(q1, 'UData', scale*qU1, 'VData', scale*qV1);
            qU2 = get(q2, 'UData');
            qV2 = get(q2, 'VData');
            set(q2, 'UData', scale*qU2, 'VData', scale*qV2);
        end
        
        function [fh] = plot_quiver_3d(obj, varargin)
            Args = {0.1};
            Args(1:(nargin-1)) = varargin;
            scale = Args{1};
            fh = figure;
            q1 = quiver3(obj.NodePositions(:,:,:,1), obj.NodePositions(:,:,:,2),...
                obj.NodePositions(:,:,:,3),...
                obj.FieldsPredicted(:,:,:,1), obj.FieldsPredicted(:,:,:,2),...
                obj.FieldsPredicted(:,:,:,3), 0);
            hold on;
            q2 = quiver3(obj.NodePositions(:,:,:,1), obj.NodePositions(:,:,:,2),...
                obj.NodePositions(:,:,:,3),...
                obj.NodeValues(:,:,:,1), obj.NodeValues(:,:,:,2), obj.NodeValues(:,:,:,3), 0);
            hold off;
            
            legend('Predicted', 'Reference');

            qU1 = get(q1, 'UData');
            qV1 = get(q1, 'VData');
            qW1 = get(q1, 'WData');
            set(q1, 'UData', scale*qU1, 'VData', scale*qV1, 'WData', scale*qW1);
            qU2 = get(q2, 'UData');
            qV2 = get(q2, 'VData');
            qW2 = get(q2, 'WData');
            set(q2, 'UData', scale*qU2, 'VData', scale*qV2, 'WData', scale*qW2);
        end
    end
    
end

