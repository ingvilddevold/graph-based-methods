classdef PartitionNet < NetworkModel
    properties
        partition;  % partition vector
        fModel;     % fine model
        fSchedule;  % fine schedule
        fState0;    % fine state0
    end

    methods
        function obj = PartitionNet(model, schedule, state0, p, fModel, fSchedule, fState0)
            obj = obj@NetworkModel(model, schedule, state0);
            obj.partition = p;

            obj.fModel = fModel;
            obj.fSchedule = fSchedule;
            obj.fState0 = fState0;
        end

        %=================================================================%
        %                            COARSENING                           %
        %=================================================================%
        function obj = mergeBlocks(obj, edges)
            % having selected some edges, 
            % merge the blocks / nodes on 
            % each side
            p = obj.partition;
            for i = 1:numel(edges)
                [s,t] = obj.network.getEndNodes(edges(i));
                if any(ismember(obj.wellNodes, s)) || any(ismember(obj.wellNodes, t))
                    warning('Selected edge contains well node and is not removed')
                    continue
                end
                p(p==t) = s;
            end
            p = processPartition(obj.fModel.G, p);
            p = compressPartition(p);            
            obj.partition = p;

            % Update model dimensions.
            obj = obj.updateFromPartition();
        end

        %------------------------------------------------------------------
        function edges = selectEdges(obj, ind, varargin)
            % Select edges across which we want to merge blocks/nodes.
            opt = struct('selectionFrac',  0.1); % coarsen top 0.1% 
            opt = merge_options(opt, varargin{:});
            edges = maxk(ind, opt.selectionFrac*size(ind));
        end

        %=================================================================%
        %                            REFINEMENT                           %
        %=================================================================%

        function obj = splitBlock(obj, block)
            p = obj.partition;
            N = obj.numBlocks();
            fineCells = find(p==block);
            for i = 1:numel(fineCells) % iterate cells in block
                newBlockNo = N + i;
                p(fineCells(i)) = newBlockNo;
            end
            p = processPartition(obj.fModel.G, p);
            p = compressPartition(p);

            %Update model dimensions.
            obj = obj.updateFromPartition(p);
        end

        %=================================================================%
        %                             UTILS                               %
        %=================================================================%
        
        function obj = updateFromPartition(obj, p)
            % Create updated model, schedule, state0 and network
            % from partition vector p

            % Model
            model = upscaleModelTPFA(obj.fModel, p);
            model = model.removeStateFunctionGroupings();
            obj.model = model;
            obj.model = obj.reimposeRelpermScaling();

            % Schedule
            obj.schedule = upscaleSchedule(obj.model, obj.fSchedule, ...
                                'wellUpscaleMethod', 'sum');
            
            % State0
            obj.state0 = upscaleState(obj.model, obj.fModel, obj.fState0);

            % Network
            obj.network = BaseGraph(model.operators.N);
            obj = obj.assignCentroidCoordsToNodes();
        end

        %------------------------------------------------------------------
        function n = numBlocks(obj)
            n = numel(unique(obj.partition));
        end
    end
end
