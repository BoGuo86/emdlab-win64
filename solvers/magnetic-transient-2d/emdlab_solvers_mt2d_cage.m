% EMDLAB: Electrical Machines Design Laboratory
% two-dimensional magnetic-transient cage
% a cage is made of a number of solid coil arms connected in parallel

classdef emdlab_solvers_mt2d_cage < handle & matlab.mixin.SetGet

    properties

        % coil arm names: all coil arms are solid
        coilArms (1,:) string;

        % resistance between adjacent bars
        Re (1,1) double {mustBeNonnegative};

        % leakage inductance between adjacent bars
        Le (1,1) double {mustBeNonnegative};

        % cage type
        type (1,:) char = 'even';

        % coils and coil arm indices
        caiStart (1,1) double;
        caiEnd (1,1) double;
        ciStart (1,1) double;
        ciEnd (1,1) double;

        % current and coltage coefficient matrices
        Ku (:,:) double;
        Kr (:,:) double;
        Kl (:,:) double;

    end

    methods

        function set.type(obj, value)

            value = lower(value);
            if ismember(value, {'even', 'odd', 'ladder'})
                obj.type = value;
            else
                error('Cage type must be <even> or <odd>');
            end

        end

        function updateKuKrKl(obj)
            n = numel(obj.coilArms);
            obj.Ku = zeros(n,n);
            for i = 1:n
                obj.Ku(i,i) = -1;
            end
            for i = 2:n
                obj.Ku(i,i-1) = 1;
            end
            switch obj.type
                case 'even'
                    obj.Ku(1,end) = 1;

                case 'odd'
                    obj.Ku(1,end) = -1;

            end
            obj.Ku = obj.Ku' *  obj.Ku;

            obj.Kr = 2 * obj.Re * eye(n);
            obj.Kl = 2 * obj.Le * eye(n);

            if (obj.Re == 0) && (obj.Le == 0)
                obj.Ku(end,:) = 0;
                obj.Kr(end,:) = 1;
            end

        end

    end

end