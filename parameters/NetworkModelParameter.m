classdef NetworkModelParameter < ModelParameter
    properties
        fillStrategy;
        mapTo;
        sens = [];
    end

    methods
        function p = NetworkModelParameter(setup, varargin)
            % Network model parameter
            % SEE ALSO: ModelParameter.m
            % setup: SimulatorSetup and NetworkModel both work
            
            opt = struct('fillStrategy', 'avg_neighbors', ...
                         'mapTo',       []);
            [opt, varargin] = merge_options(opt, varargin{:});
            
            p = p@ModelParameter(setup, varargin{:});

            p.fillStrategy = opt.fillStrategy;
            p.mapTo = opt.mapTo;
        end
    end
end