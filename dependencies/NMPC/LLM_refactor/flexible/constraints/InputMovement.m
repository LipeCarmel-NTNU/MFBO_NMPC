classdef InputMovement < handle
    % InputMovement  |Δu(k)| ≤ dumax. Rebuilt each solve (needs u_prev).
    %
    %   Emits 2*m*nu rows in physical-w coordinates; the orchestrator
    %   rescales columns to scaled-w coordinates.

    properties (SetAccess = private)
        dumax
    end

    methods
        function obj = InputMovement(dumax)
            arguments
                dumax (1,:) double {mustBePositive}
            end
            obj.dumax = dumax;
        end

        function [A, b] = linear_per_call(obj, layout, ctx)
            sz = layout.size_of('u');   % [m, nu]
            mh = sz(1);
            nu_ = sz(2);
            assert(numel(obj.dumax) == nu_, 'InputMovement:nu', ...
                'dumax size mismatch with nu.');

            cols   = layout.total_size();
            u_lin  = layout.linear_indices('u');   % m×nu
            u_prev = ctx.u_prev;

            A = zeros(2 * mh * nu_, cols);
            b = zeros(2 * mh * nu_, 1);

            row = 0;
            for k = 1 : mh
                for j = 1 : nu_
                    col_k = u_lin(k, j);
                    if k == 1
                        row = row + 1;
                        A(row, col_k) =  1;
                        b(row)        =  obj.dumax(j) + u_prev(j);

                        row = row + 1;
                        A(row, col_k) = -1;
                        b(row)        =  obj.dumax(j) - u_prev(j);
                    else
                        col_km1 = u_lin(k - 1, j);
                        row = row + 1;
                        A(row, col_k)   =  1;
                        A(row, col_km1) = -1;
                        b(row)          =  obj.dumax(j);

                        row = row + 1;
                        A(row, col_k)   = -1;
                        A(row, col_km1) =  1;
                        b(row)          =  obj.dumax(j);
                    end
                end
            end
        end
    end
end
