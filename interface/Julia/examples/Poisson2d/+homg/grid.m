classdef grid < handle
    %GRID A single grid in a multigrid heirarchy
    % compare with mgm_grid structure in mgm_multigrid.h
    
    properties
        level
        is_finest
        eig_max
        eig_min
        k_evec
        k_lam
        jacobi_omega
        jacobi_invdiag
        jacobi_inv_l1_diag
        jacobi_inv_block_diag
        gs_G
        gs_c
        ssor_M
        ssor_N
        sor_G
        sor_c
        sor_omega
        smoother
        linear_smoother
        K
        I
        Q
        L
        K_lin
        K_CG
        % Null
        % Zero
        Boundary
        Ud
        M
        R
        P
        Mesh
        Coarse  % handle to coarse grid
        debug
        dbg_spaces
        
        % variables needed for hdg solve
        refel
        is_hDG_solve
        SkelInterior2All
        SkelAll2Interior
        Bmaps
        LIFT
        VtoF
        
        pre_smooth
    end % properties
    
    methods
        function grid = grid(mesh, order, coarse)
            if ((nargin < 3) || isempty(coarse))
                grid.level = 0;
                grid.Coarse = [];
            else
                grid.level = coarse.level + 1;
                grid.Coarse = coarse;
            end
            
            grid.dbg_spaces = '      ';
            grid.dbg_spaces = grid.dbg_spaces(1:end-4*grid.level);
            grid.debug = 0;
            
            grid.Mesh = mesh;
            
            mesh.set_order(order);
            grid.sor_omega = 1;
            if (~ isempty(grid.Coarse) )
                grid.P = grid.Coarse.Mesh.assemble_interpolation(order);
                grid.R = grid.P';
            end
            
            % for Dirichlet boundary conditions
            grid.Boundary = mesh.get_boundary_node_indices(order);
            
            %% defaults ...
            grid.smoother = 'sor';
            grid.jacobi_omega = 2/3;
            grid.linear_smoother = false;
            grid.is_finest       = false;
            
            grid.is_hDG_solve = false;
        end
        
        function assemble_poisson(grid, mu)
            % fine grid material props ...
            if isnumeric(mu)
                grid.Mesh.set_muvec (mu) ;
            else
                grid.Mesh.set_coeff (mu) ;
            end
            % assemble for this level ...
            [grid.K, grid.M, grid.jacobi_inv_block_diag] = ...
                grid.Mesh.assemble_poisson(grid.Mesh.order);
            try
              syms x y z
              if ( grid.Mesh.dim == 2 )
                fx = matlabFunction(-8*pi^2*(sin(2*pi*x) * sin(2*pi*y)));
              else
                fx =matlabFunction(-12*pi^2*(sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) ));
              end
            catch
              if ( grid.Mesh.dim == 2 )
                fx = @(x,y) -8*pi^2*(sin(2*pi*x) * sin(2*pi*y));
              else
                fx = @(x,y,z) -12*pi^2*(sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) );
              end
            end
            
            grid.L = grid.Mesh.assemble_rhs(fx, grid.Mesh.order);
            grid.L(grid.Boundary) = 0;
            
            % propagate to lower grids
            if (~ isempty(grid.Coarse) )
                if isnumeric(mu)
                    % M = grid.Coarse.Mesh.assemble_mass(1);
                    % mu_coarse =   M \ ( grid.R * grid.M * mu );
                    % max(mu)
                    % max(mu_coarse)
                    harmonic = 0;
                    if (grid.Mesh.dim == 2)
                        mu2 = reshape(mu, grid.Mesh.nelems(1), grid.Mesh.nelems(2));
                        if (harmonic)
                            mu_coarse = 4 ./ ( 1./mu2(1:2:end, 1:2:end) + 1./mu2(2:2:end, 1:2:end) + 1./mu2(1:2:end, 2:2:end) + 1./mu2(2:2:end, 2:2:end) );
                        else
                            mu_coarse = 0.25*(mu2(1:2:end, 1:2:end) + mu2(2:2:end, 1:2:end) + mu2(1:2:end, 2:2:end) + mu2(2:2:end, 2:2:end));
                        end
                    else
                        mu3 = reshape(mu, grid.Mesh.nelems(1), grid.Mesh.nelems(2), grid.Mesh.nelems(3));
                        if (harmonic)
                            mu_coarse = 8 ./ ( 1./mu3(1:2:end, 1:2:end, 1:2:end) + 1./mu3(2:2:end, 1:2:end, 1:2:end) ...
                                + 1./mu3(1:2:end, 2:2:end, 1:2:end) + 1./mu3(2:2:end, 2:2:end, 1:2:end) ...
                                + 1./mu3(1:2:end, 1:2:end, 2:2:end) + 1./mu3(2:2:end, 1:2:end, 2:2:end) ...
                                + 1./mu3(1:2:end, 2:2:end, 2:2:end) + 1./mu3(2:2:end, 2:2:end, 2:2:end) );
                        else
                            mu_coarse = 0.125*(mu3(1:2:end, 1:2:end, 1:2:end) + mu3(2:2:end, 1:2:end, 1:2:end) + mu3(1:2:end, 2:2:end, 1:2:end) + mu3(2:2:end, 2:2:end, 1:2:end) + ...
                                mu3(1:2:end, 1:2:end, 2:2:end) + mu3(2:2:end, 1:2:end, 2:2:end) + mu3(1:2:end, 2:2:end, 2:2:end) + mu3(2:2:end, 2:2:end, 2:2:end) );
                        end
                    end
                    grid.Coarse.assemble_poisson (mu_coarse(:)) ;
                else
                    grid.Coarse.assemble_poisson (mu) ;
                end
            end
        end
        
        function use_linearized_smoothers(grid)
            grid.linear_smoother = true;
            [grid.K_lin, ~]  =  grid.Mesh.assemble_poisson_linearized (grid.Mesh.order);
            if (~ isempty(grid.Coarse) )
                grid.Coarse.use_linearized_smoothers();
            end
        end
        
        function set_stiffness(grid, K)
            grid.K = K;
        end
        
        % compute the residual
        function r = residual(grid, rhs, u)
            % function r = residual(grid, u, rhs)
            if ( nargin < 2 )
                rhs = grid.L;
            end
            if ( nargin < 3 )
                u = zeros(size(rhs));
            end
            
            % if ( grid.is_hDG_solve && size(u,1) ~= size(grid.K, 2))
            %    % u_hat = grid.VtoF * u;
            %    u_hat = grid.extract_skeletal_data(u);
            %    rhs_hat = grid.extract_skeletal_data(rhs);
            %    r = grid.K*u_hat - rhs_hat;
            %else
                % r = grid.K*u - rhs;
                r = rhs - grid.K*u;
            % end
        end
        
        function r = residual_lin(grid, rhs, u)
            r = grid.K_lin*u - rhs;
        end
        
        function [u, rr, iter] = solve_lin_pcg(grid, num_vcyc, smoother, smooth_steps, rhs, u)
            % disp('setting smoother');
            grid.set_smoother(smoother);
            
            % disp('computing initial residual');
            r = grid.residual(rhs, u);
            rho = zeros(size(u));
            % disp('outer v-cycle');
            rho = grid.Coarse.vcycle(smooth_steps, smooth_steps, r, rho);
            p = rho;
            disp(['Initial residual is ' num2str(norm(r))]);
            disp('------------------------------------------');
            r0 = norm(r);
            for i=1:num_vcyc
                % disp(['inner v-cycle: ' num2str(i)]);
                h = grid.K * p;
                rho_res = dot (rho, r);
                alpha = rho_res / dot ( p, h );
                u = u + alpha*p;
                r = r - alpha*h;
                
                % rho_res_prev = rho_res;
                
                disp([num2str(i, '%03d\t') ': |res| = ' num2str(norm(r),'\t%8.4e')]);
                if (norm(r)/r0 < 1e-8)
                    iter = i;
                    rr = norm(r)/r0;
                    return;
                end
                
                % precondition ..
                rho = zeros(size(u)); % needed ?
                rho = grid.Coarse.vcycle(smooth_steps, smooth_steps, r, rho);
                
                beta = dot(rho, r) / rho_res ;
                p = rho + beta*p;
            end
            disp('------------------------------------------');
            iter = num_vcyc;
            rr = norm(r)/r0;
        end
        
        
        function [u, rr, iter] = solve_pcg(grid, num_vcyc, smoother, v1, v2, rhs, u)
            % disp('setting smoother');
            grid.set_smoother(smoother);
            
            % disp('computing initial residual');
            r = grid.residual(rhs, u);
            rho = zeros(size(u));
            % disp('outer v-cycle');
            rho = grid.vcycle(v1, v2, r, rho);
            p = rho;
            disp(['Initial residual is ' num2str(norm(r))]);
            disp('------------------------------------------');
            r0 = norm(r);
            for i=1:num_vcyc
                % disp(['inner v-cycle: ' num2str(i)]);
                h = grid.K * p;
                rho_res = dot (rho, r);
                alpha = rho_res / dot ( p, h );
                u = u + alpha*p;
                r = r - alpha*h;
                
                % rho_res_prev = rho_res;
                
                disp([num2str(i, '%03d\t') ': |res| = ' num2str(norm(r),'\t%8.4e')]);
                if (norm(r)/r0 < 1e-8)
                    iter = i;
                    rr = norm(r)/r0;
                    return;
                end
                
                % precondition ..
                rho = zeros(size(u)); % needed ?
                rho = grid.vcycle(v1, v2, r, rho);
                
                beta = dot(rho, r) / rho_res ;
                p = rho + beta*p;
            end
            disp('------------------------------------------');
            iter = num_vcyc;
            rr = norm(r)/r0;
        end
        
        function [u, rr, iter] = solve_hdg_mg (grid, num_vcyc, smoother, v1, v2, rhs, u, BData)
            grid.set_smoother(smoother);
                
            num_elems = prod(grid.Mesh.nelems);
            Nfp = grid.refel.Nrp ^ (grid.refel.dim - 1);
            Nv = grid.refel.Nrp ^ (grid.refel.dim);
            
            % BData = zeros(size(grid.Bmaps));
            % 1. compute skeletal trace
            u_hat = grid.extract_skeletal_data(u);
            b = grid.hdg_residual(u_hat, rhs, BData);
            
% $$$             %%------ strongly enforce zero on the boundary----
% $$$             foo = ones(size(u_hat)); dof = length(foo);
% $$$             foo = grid.clear_skel_boundary(foo);
% $$$             index = find(foo < 0.5);
% $$$             b(index) = 0;
% $$$             K = grid.K;
% $$$             K(index,:) = 0;
% $$$             K(:,index) = 0;
% $$$             K((index - 1) * dof + index) = 1;
% $$$             grid.K = K;
% $$$             %%-----------------------------------------------
            
            u_hat_t = grid.K \ b;
            u_hat_t = grid.K \ rhs_hat;

            % 2. iterate - vcycles 
            r = grid.residual(b, u_hat);
            
% $$$             %---------- Testing (Q_0 r, Q_0 r)_0 = (r, I_1 Q_0 r)_1-----
% $$$             BData = zeros(size(grid.Bmaps));
% $$$             Q_r = grid.skel_to_cg (r, BData);
% $$$             I_Q_r = grid.cg_to_skel(Q_r);
% $$$             norm(Q_r.' * grid.M * Q_r - r.'* I_Q_r)
% $$$             keyboard
% $$$             %-----------------------------------------------------
            
            disp(['Initial residual is ' num2str(norm(r),'\t%8.4e')]);
            disp('------------------------------------------');
            r0 = norm(r);
            
            grid.plot_hdg_skel(u_hat_t);
 %           v = caxis;
            figure
            for i=1:num_vcyc
                % grid.plot_hdg_skel(u_hat-u_hat_t);
%                caxis(v);
                getframe();
                u_hat = u_hat + grid.vcycle_hdg(v1, v2, r, zeros(size(u_hat)));
                % u_hat = grid.vcycle_hdg(v1, v2, r, u_hat);
                r = grid.residual(b, u_hat);
                grid.plot_hdg_skel(r);
                temp = u_hat - u_hat_t;
                err = sqrt(temp.' * grid.K * temp);
                disp([num2str(i, '%03d\t') ': |res| = ' num2str(norm(r),'\t%8.4e') ' -- ' num2str(err,'\t%8.4e')]);

                if (norm(r)/r0 < 1e-8)
                    iter = i;
                    rr = norm(r)/r0;
                    return;
                end
            end
            disp('------------------------------------------');
            iter = num_vcyc;
            rr = norm(r)/r0;

            % grid.clear_skel_boundary(u_hat);
            % r = grid.residual(b, u_hat);
            grid.plot_hdg_skel(u_hat - u_hat_t);
            figure
            grid.plot_hdg_skel(u_hat);
            
            % 3. local solve
            lamAll = zeros(grid.Mesh.Ns_faces * Nfp,1);
            lamAll(grid.Bmaps) = BData;
            lamAll(grid.SkelInterior2All) = u_hat;
            
            u = zeros(Nv, num_elems);
            
            for e = 1:num_elems
              [u(:,e), ~, ~] = grid.localSolver(e, lamAll, rhs); 
            end

        end
        
        function [u, rr, iter] = solve(grid, num_vcyc, smoother, v1, v2, rhs, u)
            grid.set_smoother(smoother);
            
            r = grid.residual(rhs, u);
            
            disp(['Initial residual is ' num2str(norm(r))]);
            disp('------------------------------------------');
            r0 = norm(r);
            
            for i=1:num_vcyc
                u = grid.vcycle(v1, v2, rhs, u);
                r = grid.residual(rhs, u);
                disp([num2str(i, '%03d\t') ': |res| = ' num2str(norm(r),'\t%8.4e')]);
                if (norm(r)/r0 < 1e-8)
                    iter = i;
                    rr = norm(r)/r0;
                    return;
                end
            end
            disp('------------------------------------------');
            iter = num_vcyc;
            rr = norm(r)/r0;
        end
        
        function u = vcycle(grid, v1, v2, rhs, u) % <- u and not u_hat
            % function u = vcycle(grid, v1, v2, rhs, u)
            % solve system using initial guess u, given rhs
            % with v1 pre and v2 post-smoothing steps
            % disp(['CG vcycle: order ' num2str(grid.Mesh.order) ', nelems: ' num2str(grid.Mesh.nelems(1)) 'X' num2str(grid.Mesh.nelems(2))]);
            
            if ( isempty( grid.Coarse ) )
                if (grid.linear_smoother)
                    u = grid.K_lin \ rhs;
                else
                    u = grid.K \ rhs;
%                    u = grid.K \ (grid.M * rhs);
                end
                
                return;
            end
            
            % 1. pre-smooth
            u = grid.smooth ( v1, rhs, u );
                
            % 2. compute residual
            if (grid.linear_smoother && ~grid.is_finest)
                % disp('linear residual');
                res = grid.residual_lin(rhs, u);
            else
                % disp('high-order residual');
                res = grid.residual(rhs, u);
            end
                
           % 3. restrict
           res_coarse = grid.R * res;
           res_coarse(grid.Coarse.Boundary) = 0;
            
           % 4. ---------- recurse -----------
           u_corr_coarse = grid.Coarse.vcycle(v1, v2, res_coarse, zeros(size(res_coarse)));
            
           % 5. prolong and correct
           u = u + grid.P * u_corr_coarse;
           
           % 6. post-smooth
           u = grid.smooth ( v2, rhs, u );
           % grid.plot_spectrum(u, 'g', rhs);
            
       end % v-cycle
       
       % hDG v-cycle
        function u_hat = vcycle_hdg (grid, v1, v2, rhs_hat, u_hat) % <- u_hat and not u
            % disp(['hdg vcycle: order ' num2str(grid.refel.N) ', nelems: ' num2str(grid.Mesh.nelems(1)) 'X' num2str(grid.Mesh.nelems(2))]);
            % SMOOTH
            grid.pre_smooth = 1;
            u_hat = grid.smooth ( v1, rhs_hat, zeros(size(u_hat)) );
%            u_hat = grid.clear_skel_boundary(u_hat);
            
            % Compute RESIDUAL
            res = grid.residual(rhs_hat, u_hat);
%            res = grid.clear_skel_boundary(res);

% $$$             grid.plot_hdg_skel(u_hat);
% $$$             getframe();
            
            if ( ~ grid.Coarse.is_hDG_solve )
              BData = zeros(size(grid.Bmaps));
              
              % RESTRICT
              res_coarse = grid.skel_to_cg (res, BData);
              %res_coarse(grid.Coarse.Boundary) = 0;
              
              % VCYCLE
              u_corr_coarse = grid.Coarse.vcycle(v1, v2, res_coarse, zeros(size(res_coarse)));
              
              % PROLONG + correct 
              u_hat_corr = grid.cg_to_skel (u_corr_coarse);
%              u_hat_corr = grid.clear_skel_boundary(u_hat_corr);
              u_hat = u_hat + u_hat_corr; 

            else
              % RESTRICT
              res_coarse = grid.hdg_restrict(res);

              % VCYCLE
              u_corr_coarse = grid.Coarse.vcycle_hdg(v1, v2, res_coarse, zeros(size(res_coarse)));

              % PROLONG + correct
              u_hat = u_hat + grid.hdg_prolong(u_corr_coarse);
            end

            res = grid.residual(rhs_hat, u_hat);
            
            % SMOOTH 
% $$$             u_hat = u_hat + grid.smoother_chebyshev_adjoint (v2, res, zeros(size(u_hat)) );
            u_hat = u_hat + grid.smooth( v2, res, zeros(size(u_hat)) );
%            u_hat = grid.clear_skel_boundary(u_hat);
  
% $$$             u_hat = grid.smooth ( v2, rhs_hat, u_hat );
% $$$ 
            u_hat = u_hat + grid.smoother_chebyshev_adjoint (v2, rhs_hat, u_hat );
 
        end % v-cycle
        
        % smoothers
        function u = smooth (grid, v, rhs, u)
            switch(grid.smoother)
                case 'jacobi',
                    u = grid.smoother_jacobi(v, rhs, u);
                    return;
                case 'l1_jac',
                    u = grid.smoother_l1_jacobi(v, rhs, u);
                    return;
                case 'blk_jac',
                    u = grid.smoother_block_jacobi(v, rhs, u);
                    return;
                case 'gs',
                    grid.sor_omega = 1.0;
                    u = grid.smoother_gauss_seidel(v, rhs, u);
                    return;
                case 'chebssor'
                    u = grid.smoother_chebyshev_ssor (v, rhs, u);
                    return;
                case 'chebyshev2'
                    u = grid.smoother_chebyshev2 (v, rhs, u);
                    return;
                case 'chebyshev',
                    u = grid.smoother_chebyshev(v, rhs, u);
                    return;
                case 'sor',
                    u = grid.smoother_sor(v, rhs, u);
                    return;
                case 'ssor',
                    u = grid.smoother_sym_sor(v, rhs, u);
                    return;
                case '2sr',
                    u = grid.smoother_2sr(v, rhs, u);
                    return;
                case 'hybrid',
                    u = grid.smoother_hybrid(v, rhs, u);
                    return;
                otherwise
                    disp('ERROR: Unrecognized smoother type');
                    return;
            end
        end
        
        function set_coeff(grid, mu)
            grid.Mesh.set_coeff (mu) ;
            if (~ isempty(grid.Coarse) )
                grid.Coarse.Mesh.set_coeff (mu) ;
            end
        end
        
        function set_smoother(grid, sm)
            grid.smoother = sm;
            if (~ isempty(grid.Coarse) )
                grid.Coarse.set_smoother(sm);
            end
        end
        
        function u = smoother_jacobi (grid, v, rhs, u)
            % standard jacobi smoother
            if ( isempty(grid.jacobi_invdiag) )
                % Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
                if (grid.linear_smoother && ~grid.is_finest)
                    D = diag(grid.K_lin);
                else
                    D = diag(grid.K);
                end
                grid.jacobi_invdiag = 1./D;
            end
            
            for i=1:v
                if (grid.linear_smoother)
                    res  = grid.jacobi_invdiag .* grid.residual_lin(rhs, u);
                else
                    res  = grid.jacobi_invdiag .* grid.residual(rhs, u);
                end
                % res  = grid.jacobi_invdiag .* grid.residual(rhs, u);
                u = u + grid.jacobi_omega.*res;
                % r = norm(res);
                % disp([grid.dbg_spaces num2str(r)]);
                % norm(r)
            end
        end % jacobi
        
        function u = smoother_l1_jacobi (grid, v, rhs, u)
            % l1 jacobi smoother
            if ( isempty(grid.jacobi_inv_l1_diag) )
                
                D = diag(grid.K) + sum(abs(grid.K - diag(diag(grid.K)) ), 2) ;
                
                grid.jacobi_inv_l1_diag = 1./D;
            end
            
            for i=1:v
                res  = grid.jacobi_inv_l1_diag .* grid.residual(rhs, u);
                % res  = grid.jacobi_invdiag .* grid.residual(rhs, u);
                u = u + grid.jacobi_omega.*res;
                % r = norm(res);
                % disp([grid.dbg_spaces num2str(r)]);
                % norm(r)
            end
        end % l1-jacobi
        
        function u = smoother_block_jacobi (grid, v, rhs, u)
            % block jacobi smoother
            if ( isempty(grid.jacobi_inv_block_diag) )
                error('inv block doagonal not assembled');
                % grid.jacobi_inv_block_diag = grid.assemble_inv_block();
            end
            
            for i=1:v
                if (grid.linear_smoother)
                    res  = grid.jacobi_inv_block_diag * grid.residual_lin(rhs, u);
                else
                    res  = grid.jacobi_inv_block_diag * grid.residual(rhs, u);
                end
                % res  = grid.jacobi_invdiag .* grid.residual(rhs, u);
                u = u + grid.jacobi_omega.*res;
                % r = norm(res);
                % disp([grid.dbg_spaces num2str(r)]);
                % norm(r)
            end
        end % jacobi
        
        function u = smoother_sor (grid, v, rhs, u)
            if ( isempty ( grid.sor_G ) )
                % Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
                kL = tril(grid.K, -1);
                kD = diag(diag(K));
                kU = triu(Kc, 1);
                grid.sor_G = - (kD + grid.sor_omega*kL) \ (grid.sor_omega*kU + (grid.sor_omega - 1.0)*kD );
                grid.sor_c = grid.ZeroBoundary *( (kD + grid.sor_omega*kL) \ (grid.sor_omega*rhs) );
            end
            
            for i=1:v
                u = grid.sor_G*u + grid.sor_c;
            end
        end
        
        function u = smoother_sym_sor (grid, v, rhs, u)
            if ( isempty ( grid.ssor_M ) )
                w = grid.sor_omega;
                n = length(u);
                if ( grid.linear_smoother )
                    grid.ssor_M = spdiags( (1/w)*diag(grid.K_lin), 0, n, n) + tril(grid.K_lin,-1);
                    grid.ssor_N = spdiags(((1-w)/w)*diag(grid.K_lin), 0, n, n) - triu(grid.K_lin,1);
                else
                    grid.ssor_M = spdiags( (1/w)*diag(grid.K), 0, n, n) + tril(grid.K,-1);
                    grid.ssor_N = spdiags(((1-w)/w)*diag(grid.K), 0, n, n) - triu(grid.K,1);
                end
            end
            
            for i=1:v
                if (grid.linear_smoother)
                    r = grid.residual_lin(rhs, u);
                else
                    r = grid.residual(rhs, u);
                end
                
%                 if ( grid.is_hDG_solve) 
%                    if (grid.pre_smooth)
%                        u = u + grid.ssor_M \ r;
%                    else
%                        u = grid.ssor_M' \ (grid.ssor_N'*u + rhs);
%                    end
%                 else
                    u = u + grid.ssor_M \ r;
                    u = grid.ssor_M' \ (grid.ssor_N'*u + rhs);
                    
                    res = grid.residual ( rhs, u );
                    r = norm(res);
                    disp([' ---- ' num2str(i) ' : ' num2str(r)]);
%                 end
            end
        end
        
        function u = smoother_chebyshev_jacobi (grid, v, rhs, u)
            if ( isempty ( grid.eig_max ) )
                % Kc = grid.Null' * grid.K * grid.Null;
                Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
                D = diag(Kc);
                grid.jacobi_invdiag = 1./D;
                Kc = (eye(size(Kc)) - (diag(D) \ Kc) );
                grid.eig_max = eigs(Kc, 1, 'LM');
                grid.eig_min = eigs(Kc, 1, 'SM');
            end
            
            l_max = grid.eig_max;
            l_min =  (grid.eig_min + grid.eig_max)/2;
            
            rho = 2/(l_min + l_max);
            
            mu_0 = 1;
            mu_1 = rho;
            y_0 = u;
            y_1 = grid.smoother_jacobi(1, rhs, u);
            for i=2:v
                mu_2 = 1.0 / ( 2.0/(rho*mu_1) - 1.0/mu_0);
                u = (2.0*mu_2)/(rho*mu_1)* grid.smoother_jacobi(1, rhs, y_1) - (mu_2/mu_1)*y_0;
                y_0 = y_1; mu_0 = mu_1;
                y_1 =  u ; mu_1 = mu_2;
            end
        end
        
        function u = smoother_chebyshev_ssor (grid, v, rhs, u)
            if ( isempty ( grid.eig_max ) )
                % Kc = grid.Null' * grid.K * grid.Null;
                Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
                grid.eig_max = eigs(Kc, 1, 'LM');
                grid.eig_min = eigs(Kc, 1, 'SM');
            end
            if ( isempty ( grid.ssor_M ) )
                w = grid.sor_omega;
                n = length(u);
                grid.ssor_M = spdiags( (1/w)*diag(grid.K), 0, n, n) + tril(grid.K,-1);
                grid.ssor_N = spdiags(((1-w)/w)*diag(grid.K), 0, n, n) - triu(grid.K,1);
            end
            
            l_max = grid.eig_max;
            l_min =  (grid.eig_min + grid.eig_max)/2;
            
            rho = 2/(l_min + l_max);
            
            mu_0 = 1;
            mu_1 = rho;
            y_0 = u;
            y_1 = grid.smoother_sym_sor(1, rhs, u);
            for i=2:v
                mu_2 = 1.0 / ( 2.0/(rho*mu_1) - 1.0/mu_0);
                u = (2.0*mu_2)/(rho*mu_1)* grid.smoother_sym_sor(1, rhs, y_1) - (mu_2/mu_1)*y_0;
                y_0 = y_1; mu_0 = mu_1;
                y_1 =  u ; mu_1 = mu_2;
            end
            
        end
        
        function set_sor_omega(grid, w)
            grid.sor_omega = w;
            grid.sor_G = [];
            grid.sor_c = [];
            grid.ssor_M = [];
            grid.ssor_N = [];
        end
        
        function u = smoother_gauss_seidel (grid, v, rhs, u)
            if ( isempty ( grid.gs_G ) )
                Kc = grid.K; %(eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
                LD = tril(Kc);
                grid.gs_G = -LD \ triu(Kc, 1);
                grid.gs_c = ( LD \ rhs );
                grid.gs_c(grid.Boundary) = 0;
            end
            
            for i=1:v
                if (grid.pre_smooth)
                    u = grid.gs_G*u + grid.gs_c;
                else
                    u = grid.gs_G' * u + grid.gs_c;
                end
            end
        end
        
        function u = smoother_gauss_seidel_adjoint (grid, v, rhs, u)
            if ( isempty ( grid.gs_G ) )
                Kc = fliplr(grid.K); %(eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
                LD = tril(Kc);
                grid.gs_G = -LD \ triu(Kc, 1);
                grid.gs_c = ( LD \ rhs );
                grid.gs_c(grid.Boundary) = 0;
            end
            
            for i=1:v
                u = grid.gs_G*u + grid.gs_c;
            end
            u = flipud(u);
        end
        
%         function u = smoother_gauss_seidel_adjoint (grid, v, rhs, u)
%             if ( isempty ( grid.gs_G ) )
%                 Kc = fliplr(grid.K); %(eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
%                 LD = tril(Kc);
%                 grid.gs_G = -LD \ triu(Kc, 1);
%                 grid.gs_c = ( LD \ rhs );
%                 grid.gs_c(grid.Boundary) = 0;
%             end
%             
%             for i=1:v
%                 u = grid.gs_G*u + grid.gs_c;
%             end
%             u = flipud(u);
%         end
        
        function u = smoother_hybrid (grid, v, rhs, u)
            % u = grid.smoother_gauss_seidel (v, rhs, u);
            u = grid.smoother_chebyshev (v, rhs, u);
            % u = grid.smoother_chebyshev (v, rhs, u);
            u = grid.smoother_sym_sor(2, rhs, u);
        end
        
        function u = smoother_2sr (grid, v, rhs, u)
            % 2-step stationary iterative smoother
            % factors
            if ( isempty ( grid.eig_max ) )
                % Kc = grid.Null' * grid.K * grid.Null;
                Kc = (eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
                grid.eig_max = eigs(Kc,1, 'LM');
                grid.eig_min = eigs(Kc,1, 'SM');
            end
            
            l_max = grid.eig_max;
            % l_min = grid.eig_max*.9;
            l_min = (grid.eig_min + grid.eig_max)/2;
            
            rho       = (1 - l_min/l_max)/(1 + l_min/l_max);
            alpha     = 2/( 1 + sqrt(1-rho*rho));
            epsilon   = 2/(l_min + l_max);
            epsalpha  = epsilon * alpha;
            
            % variables
            u0  = zeros(size(u));
            res = -rhs; % grid.residual(rhs, u);
            u1 = epsilon * res ;
            
            for iter = 1:v,                            % begin iteration
                res = grid.residual ( rhs, u );
                u = alpha*u1 + (1-alpha)*u0 - epsalpha*res;
                u0 = u1; u1 = u;
                % r = norm(res);
                % n0 = norm(u0);
                % n1 = norm(u1);
                % disp([grid.dbg_spaces 'residual: ' num2str(r)]); % ' u0: ' num2str(n0) ' u1: ' num2str(n1)]);
            end % end iteration
            
        end % 2sr
        
        function u = smoother_chebyshev2 (grid, v, rhs, u)
            if ( isempty ( grid.eig_max ) )
                disp('computing eigenvalues');
                tic;
                % Kc = grid.Null' * grid.K * grid.Null;
                Kc = grid.K; %(eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
                % d = eigs(Kc, 2, 'be');
                grid.eig_max = eigs(Kc, 1, 'lm');
                grid.eig_min = eigs(Kc, 1, 'sm');
                toc;
            end
            
            % adjust the eigenvalues to hit the upper spectrum
            l_max = grid.eig_max;
            l_min = (grid.eig_min + grid.eig_max)/2;
            
            c = (l_min - l_max)/2;
            d = (l_min + l_max)/2;
            
            p = zeros(size(u));
            
            for iter = 1:v
                res = grid.residual ( rhs, u );
                %r = norm(res)
                % disp([grid.dbg_spaces 'residual: ' num2str(r)]);
                if ( iter == 1 )
                    alpha = 1.0/d;
                elseif (iter == 2)
                    alpha = 2*d / (2*d*d - c*c);
                else
                    alpha = 1.0/(d - alpha*c*c*0.25);
                end
                
                beta = alpha * d - 1.0;
                
                p = -alpha * res + beta * p;
                u = u + p;
            end
        end % chebyshev
        
        function u = smoother_chebyshev (grid, v, rhs, u)
            if ( isempty ( grid.eig_max ) )
                if ( grid.linear_smoother )
                    D = diag(grid.K_lin);
                    grid.jacobi_invdiag = 1./D;
                    Kc = spdiags(grid.jacobi_invdiag,0,length(D), length(D)) * grid.K;
                else
                    D = diag(grid.K);
                    grid.jacobi_invdiag = 1./D;
                    Kc = spdiags(grid.jacobi_invdiag,0,length(D), length(D)) * grid.K;
                end
                opts.tol = 0.01;
                grid.eig_max = eigs(Kc, 1, 'lm', opts);
                % grid.eig_min = eigs(Kc, 1, 'sm');
            end
            
            % adjust the eigenvalues to hit the upper spectrum
            beta = grid.eig_max;
            alpha = 0.25*grid.eig_max;% (grid.eig_min + grid.eig_max)/2;
            
            delta = (beta - alpha)/2;
            theta = (beta + alpha)/2;
            s1 = theta/delta;
            rhok = 1./s1;
            
            d = zeros(size(u));
            
            % first loop
            if (grid.linear_smoother && ~grid.is_finest)
                res = grid.residual_lin ( rhs, u );
            else
                res = grid.residual ( rhs, u );
            end
            d = res/theta.* grid.jacobi_invdiag;
            u = u + d;
            
            for iter = 2:v
                rhokp1 = 1/ (2*s1 - rhok);
                d1 = rhokp1 * rhok;
                d2 = 2*rhokp1 / delta;
                rhok = rhokp1;
                if (grid.linear_smoother && ~grid.is_finest)
                    res = grid.residual_lin ( rhs, u );
                else
                    res = grid.residual ( rhs, u );
                end
                disp(['  -- ' num2str(iter) ' : ' num2str(norm(res))]);
                d = d1 * d + d2 * res.*grid.jacobi_invdiag;
                u = u + d;
            end
        end % chebyshev

        function u = smoother_chebyshev_adjoint (grid, v, rhs, u)
            if ( isempty ( grid.eig_max ) )
              K = fliplr(grid.K);
              D = diag(K);
              grid.jacobi_invdiag = 1./D;
              Kc = spdiags(grid.jacobi_invdiag,0,length(D), length(D)) * K;
              
              opts.tol = 0.01;
              grid.eig_max = eigs(Kc, 1, 'lm', opts);
                % grid.eig_min = eigs(Kc, 1, 'sm');
            end
            
            % adjust the eigenvalues to hit the upper spectrum
            beta = grid.eig_max;
            alpha = 0.25*grid.eig_max;% (grid.eig_min + grid.eig_max)/2;
            
            delta = (beta - alpha)/2;
            theta = (beta + alpha)/2;
            s1 = theta/delta;
            rhok = 1./s1;
            
            d = zeros(size(u));
            
            % first loop
            if (grid.linear_smoother && ~grid.is_finest)
                res = grid.residual_lin ( rhs, u );
            else
                res = grid.residual ( rhs, u );
            end
            d = res/theta.* grid.jacobi_invdiag;
            u = u + d;
            
            for iter = 2:v
                rhokp1 = 1/ (2*s1 - rhok);
                d1 = rhokp1 * rhok;
                d2 = 2*rhokp1 / delta;
                rhok = rhokp1;
                if (grid.linear_smoother && ~grid.is_finest)
                    res = grid.residual_lin ( rhs, u );
                else
                    res = grid.residual ( rhs, u );
                end
                % disp([num2str(iter) ':' num2str(norm(res))]);
                d = d1 * d + d2 * res.*grid.jacobi_invdiag;
                u = u + d;
            end

            u = flipud(u);
        end % chebyshev

%         function u = smoother_chebyshev_adjoint (grid, v, rhs, u)
%             if ( isempty ( grid.eig_max ) )
%               K = fliplr(grid.K);
%               D = diag(K);
%               grid.jacobi_invdiag = 1./D;
%               Kc = spdiags(grid.jacobi_invdiag,0,length(D), length(D)) * K;
%               
%               opts.tol = 0.01;
%               grid.eig_max = eigs(Kc, 1, 'lm', opts);
%                 % grid.eig_min = eigs(Kc, 1, 'sm');
%             end
%             
%             % adjust the eigenvalues to hit the upper spectrum
%             beta = grid.eig_max;
%             alpha = 0.25*grid.eig_max;% (grid.eig_min + grid.eig_max)/2;
%             
%             delta = (beta - alpha)/2;
%             theta = (beta + alpha)/2;
%             s1 = theta/delta;
%             rhok = 1./s1;
%             
%             d = zeros(size(u));
%             
%             % first loop
%             if (grid.linear_smoother && ~grid.is_finest)
%                 res = grid.residual_lin ( rhs, u );
%             else
%                 res = grid.residual ( rhs, u );
%             end
%             d = res/theta.* grid.jacobi_invdiag;
%             u = u + d;
%             
%             for iter = 2:v
%                 rhokp1 = 1/ (2*s1 - rhok);
%                 d1 = rhokp1 * rhok;
%                 d2 = 2*rhokp1 / delta;
%                 rhok = rhokp1;
%                 if (grid.linear_smoother && ~grid.is_finest)
%                     res = grid.residual_lin ( rhs, u );
%                 else
%                     res = grid.residual ( rhs, u );
%                 end
%                 % disp([num2str(iter) ':' num2str(norm(res))]);
%                 d = d1 * d + d2 * res.*grid.jacobi_invdiag;
%                 u = u + d;
%             end
% 
%             u = flipud(u);
%         end % chebyshev

        
        function [evec, eval] = get_eigenvectors(grid)
            % generate the correct matrix
            Kc = grid.K; %(eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
            [evec, eval] = svd(full(Kc)); %eig(full(Kc), full(grid.M));
            [eval,per] = sort(diag(eval),'ascend');
            evec = evec(:,per);
            grid.k_evec = evec;
            grid.k_lam = eval;
        end
        
        function plot_spectrum(g, u, clr, rhs)
            if (g.debug)
%                 subplot(1,2,1);
                a = u; % g.M * 
                q = repmat(a, 1, 80);
                b = abs(dot (g.k_evec, q));
                % plot eigenvalues
                plot(b, clr); hold on;
%                 subplot(1,2,2);
%                 rr = g.residual(rhs, u);
%                 n = sqrt(length(rr));
%                 imagesc(reshape(rr, n, n)); colorbar; hold off;
%                 % grid on;
%                 odr = g.Mesh.fem.shape;
%                 set(gca, 'xtick', odr+0.5:odr:odr*g.Mesh.nelem);
%                 set(gca, 'ytick', odr+0.5:odr:odr*g.Mesh.nelem);
            end
        end
        
        function u0 = get_u0(grid)
            if (grid.debug)
                if ( isempty( grid.k_evec ) )
                    [grid.k_evec, ~] = eigs(grid.K, grid.M, 80, 'BE');
                end
                n = size(grid.k_evec, 2);
                lam = ones(n,1);
                % lam(1:n/4) = 1;
                u0 = grid.k_evec*lam;
            else
                u0 = rand(size(grid.L()));
                u0(grid.Boundary) = 0;
            end
        end
        
        %%~~~~~~~~~~ hDG functions ~~~~~~~~~~
        %
        % to use multigrid, use solve or solve_pcg instead
        function [u, qx, qy] = solve_hdg (self, forcing, Bdata)
            % full solve without multigrid,
            disp('in solve HDG');
            K = prod(self.Mesh.nelems);
            Nfp = self.refel.Nrp ^ (self.refel.dim - 1);
            Nv = self.refel.Nrp ^ (self.refel.dim);
            
            lamInterior = zeros(size(self.SkelInterior2All));
            
            rhs = self.hdg_residual(lamInterior, forcing, Bdata);
            
            lamInterior = self.K \ rhs;
            
            lamAll = zeros(self.Mesh.Ns_faces * Nfp,1);
            lamAll(self.Bmaps) = Bdata;
            lamAll(self.SkelInterior2All) = lamInterior;
            
            u = zeros(Nv,K);
            qx = zeros(Nv,K);
            qy = zeros(Nv,K);
            
            for e = 1:K
                [u(:,e),qx(:,e),qy(:,e)] = self.localSolver(e, lamAll, forcing);
            end
            
            disp('leaving solve HDG');
        end
        
        function u_hat = extract_skeletal_data(self, u)
            u_hat = zeros(size(self.SkelInterior2All));
            
            % first approximation - average ...
            for  sf=1:self.Mesh.Ns_faces
                [e1, f1, e2, f2]  = self.Mesh.get_face_elements(sf);
                
                if (e1 > 0) && (e2 > 0) % interior faces
                    if (f1 > 0)
                        pts = self.Mesh.element_nodes(e1, self.refel);
                        idxf = self.Mesh.get_skeletal_face_indices(self.refel, e1, f1);
                        idxv = self.Mesh.get_discontinuous_face_indices(self.refel, e1, f1);
                        u_hat(self.SkelAll2Interior(idxf)) = u_hat(self.SkelAll2Interior(idxf)) + u(idxv);
                    end
                    if (f2 > 0)
                        pts = self.Mesh.element_nodes(e2, self.refel);
                        idxf = self.Mesh.get_skeletal_face_indices(self.refel, e2, f2);
                        idxv = self.Mesh.get_discontinuous_face_indices(self.refel, e2, f2);
                        u_hat(self.SkelAll2Interior(idxf)) = u_hat(self.SkelAll2Interior(idxf)) + u(idxv);
                    end
                end
            end
            u_hat = 0.5 .* u_hat;
        end
        
        function [u_skel, Bdata] = cg_to_skel (self, u_cg)
            Nfp = self.refel.Nrp ^ (self.refel.dim - 1);
            num_elem = prod(self.Mesh.nelems);
            
            lamAll = zeros(self.Mesh.Ns_faces * Nfp,1);
            
            % loop over elements
            for e=1:num_elem
                % TODO: remove hard-coding with 4 here
                for fid=1:4
                    idx = self.Mesh.get_skeletal_face_indices(self.refel, e, fid);
                    [cg_idx, gfid] = self.Mesh.get_continuous_face_indices(self.refel, e, fid);
                    
                    % Prolongate CG linear solution to HDG linear skeleton
                    lamAll(idx) = u_cg(cg_idx);
                end
            end

            % lamAll(self.Bmaps) = 0;
            
            Bdata  = lamAll(self.Bmaps);
            u_skel = lamAll(self.SkelInterior2All);
            
% $$$             Bdata  = zeros(length(self.Bmaps),1);
% $$$             u_skel = self.I'*u_cg;
        end
                
        function u_cg = skel_to_cg (self, u_skel, Bdata)

            Nfp = self.refel.Nrp ^ (self.refel.dim - 1);
            
            dof = prod(self.Mesh.nelems*self.refel.N + 1);
            
            num_elem = prod(self.Mesh.nelems);
            
            nx = self.Mesh.nelems(1);
            ny = self.Mesh.nelems(2);
            
            % TODO: need to replace 1.0 with the actual domain length
            hx = 1.0/nx; hy = 1.0/ny;
            fac = hx*hy / (2*(hx+hy));
            
            % add Bdata to skel to get linear dg
            lamAll = zeros(self.Mesh.Ns_faces * Nfp,1);
            lamAll(self.Bmaps) = Bdata;
            lamAll(self.SkelInterior2All) = u_skel;
            
            u_cg = zeros(dof,1);
            
            
            % loop over elements
            for e=1:num_elem
                % TODO: remove hard-coding with 4 here
                for fid=1:4
                    idx = self.Mesh.get_skeletal_face_indices(self.refel, e, fid);
                    [cg_idx, gfid] = self.Mesh.get_continuous_face_indices(self.refel, e, fid);
                    
                    %----Restrict HDG linear skeleton solution to linear CG mesh-------
                    % (Q_0 mu, v)_0 = (mu, I_1 v)_1
% $$$                     Jf = self.Mesh.geometric_factors_face(self.refel, e, fid);
% $$$                     Md = self.refel.w .* Jf;
% $$$                     Mf = self.refel.q1d' * diag(Md) * self.refel.q1d;
% $$$                     rhs = Mf * lamAll(idx);  
                    rhs = lamAll(idx);
                    % rhs = Jf .* (self.refel.Mr * lamAll(idx));

                    if self.SkelAll2Interior(idx(1)) < eps
                      if norm(rhs) > 1.e-14
                        error('rhs must be zero here');
                      end
                        u_cg(cg_idx) = u_cg(cg_idx) + rhs;
                    else
                        u_cg(cg_idx) = u_cg(cg_idx) + 0.5 * rhs;
                    end
                    %------------------------------------------------------------------
                end
            end
            % u_cg = u_cg';
            if ( isempty(self.M) )
                self.M = self.Mesh.assemble_mass(self.refel.N);
            end
% $$$             
% $$$             bdy_index = ...
% $$$                 self.Mesh.get_boundary_node_indices(self.Mesh.order);
% $$$             u_cg(bdy_index)  = 0; 
% $$$         
% $$$ %                u_cg = self.I * u_skel;
% $$$ 
% $$$ 
% $$$             M = self.M;
% $$$             
% $$$             M(bdy_index,:) = 0;
% $$$             M(:, bdy_index) = 0;
% $$$             M((bdy_index - 1) * dof + bdy_index) = 1;
% $$$             
% $$$             self.M = M;
            
%            u_cg = self.M \ u_cg;
        end
        
        
        % restrict from fine to coarse hDG skel ... 
        function uhat_coarse = hdg_restrict(self, uhat_fine)
            Nrp_fine     = self.refel.Nrp;
            Nrp_coarse   = self.Coarse.refel.Nrp;
            
            num_faces = length(uhat_fine) / Nrp_fine;
            
            uhat_coarse = zeros(num_faces * Nrp_coarse, 1);
            
            p = self.Coarse.refel.p_p_1d;
            r = p';
 
            for f=1:num_faces
              % restrict from fine to coarse face
              uhat_coarse ( ((f-1)*Nrp_coarse+1):(f*Nrp_coarse) ) = r * uhat_fine ( ((f-1)*Nrp_fine+1):(f*Nrp_fine) );
            end
        end
        
        function uhat_fine = hdg_prolong(self, uhat_coarse)
            Nrp_fine     = self.refel.Nrp;
            Nrp_coarse   = self.Coarse.refel.Nrp;
            
            num_faces = length(uhat_coarse) / Nrp_coarse;
            
            uhat_fine = zeros(num_faces * Nrp_fine, 1);
            
            p = self.Coarse.refel.p_p_1d;
 
            for f=1:num_faces
              % prolong from coarse to fine face
              uhat_fine ( ((f-1)*Nrp_fine+1):(f*Nrp_fine) ) = p * uhat_coarse ( ((f-1)*Nrp_coarse+1):(f*Nrp_coarse) );
            end
        end
        
        function res = hdg_residual(self, lamInterior, forcing, Bdata)
            % InteriorF2AllF = HDG.InteriorF2AllF;
            
            refel = self.refel;
            K = prod(self.Mesh.nelems);
            Nfp = self.refel.Nrp ^ (self.refel.dim - 1);
            Nfaces = self.refel.dim * 2;
            
            nx   = self.Mesh.nx;
            ny   = self.Mesh.ny;
            taur = self.Mesh.taur;
            
            res = zeros(Nfp * self.Mesh.Ni_faces, 1);
            
            % lamAll
            lamAll = zeros(self.Mesh.Ns_faces * Nfp, 1);
            lamAll(self.Bmaps) = Bdata;
            lamAll(self.SkelInterior2All) = lamInterior;
            
            for e = 1:K
                % e1 solution
                [u, qx, qy] = self.localSolver(e, lamAll, forcing);
                
                for f = 1:Nfaces
                    idxf = self.Mesh.get_skeletal_face_indices(refel, e, f);
                    
                    if abs( self.SkelAll2Interior( idxf(1) ) ) < eps
                        continue;
                    end
                    lam = lamAll(idxf);
                    
                    Jf   = self.Mesh.geometric_factors_face(refel,e,f);
                    idxv = self.Mesh.get_discontinuous_face_indices(refel, 1, f);
                    
                    
                    % construct the residual
                    Fhat = Jf .* (refel.Mr * (qx(idxv) * nx(f) + qy(idxv) * ny(f) + ...
                        taur * (u(idxv) - lam)));
                    
                    res(self.SkelAll2Interior(idxf)) = res(self.SkelAll2Interior(idxf)) + Fhat;
                end
            end
        end
        
        function [du_dlam, dqx_dlam, dqy_dlam] = DifflocalSolver(self, e1)
            % compute derivatives for the local solver
            
            % predefined normal vector, don't like it but stick with it for now
            nx   = self.Mesh.nx;
            ny   = self.Mesh.ny;
            taur = self.Mesh.taur; % maybe move to grid
            
            Nfp = self.refel.Nrp ^ (self.refel.dim - 1);
            Nv = self.refel.Nrp ^ (self.refel.dim);
            
            uu  = zeros(Nv,Nv);
            uqx = zeros(Nv,Nv);
            uqy = zeros(Nv,Nv);
            
            pts = self.Mesh.element_nodes(e1, self.refel);
            [Jv, Dv] = self.Mesh.geometric_factors(self.refel, pts);
            
            eMat = self.Mesh.element_mass(e1, self.refel, Jv);
            eMatInv = inv(eMat);
            
            % number of faces
            Nfaces = self.refel.dim * 2;
            
            du_dlam  = zeros(Nv, Nfp * Nfaces);
            dqx_dlam = zeros(Nv, Nfp * Nfaces);
            dqy_dlam = zeros(Nv, Nfp * Nfaces);
            
            drhsqx_dlam = zeros(Nv, Nfp * Nfaces);
            drhsqy_dlam = zeros(Nv, Nfp * Nfaces);
            drhu_dlam   = zeros(Nv, Nfp * Nfaces);
            
            dF_dlam     = zeros(Nv, Nfp * Nfaces);
            
            % advection stiffness
            [Kex, Key] = self.Mesh.element_stiffness_advection (e1, self.refel, Jv, Dv);
            
            uqx = -Kex;
            uqy = -Key;
            
            % residual for qx and qy equations
            for f = 1:Nfaces %
                index = (f-1)*Nfp+1:f*Nfp;
                
                % geometrix factors at gll points on face
                Jf = self.Mesh.geometric_factors_face(self.refel,e1,f);
                
                %-- residual due to lambda
                % rhsfx = Jf .* (refel.Mr * lamlocal) * nx(f);
                %-- lift to volume residual q equation
                % rhsqx(idxv) = rhsqx(idxv) + rhsfx;
                drhsqx_dlam(:,index) = -self.LIFT(:,:,f) * (diag(Jf .* nx(f)));
                
                % rhsfy = Jf .* (refel.Mr * lamlocal) * ny(f);
                % rhsqy(idxv) = rhsqy(idxv) + rhsfy;
                drhsqy_dlam(:,index) = -self.LIFT(:,:,f) * (diag(Jf .* ny(f)));
                
                % lift to volume residual u equation
                % rhsu(idxv)  = rhsu(idxv) - ...
                %    taur * Jf .* (refel.Mr * lamlocal);
                drhsu_dlam(:,index) = self.LIFT(:,:,f) * (diag(Jf .* taur));
                
                % lift to volume for uu
                bdry =  self.LIFT(:,:,f) * (diag(Jf) * self.VtoF(:,:,f));
                bdryx = self.LIFT(:,:,f) * (diag(Jf) * self.VtoF(:,:,f)) * nx(f);
                bdryy = self.LIFT(:,:,f) * (diag(Jf) * self.VtoF(:,:,f)) * ny(f);
                
                uu  = uu  + taur * bdry;
                uqx = uqx +        bdryx;
                uqy = uqy +        bdryy;
            end
            
            qxMatrix = eMatInv * Kex;
            % qxrhs = eMatInv * rhsqx;
            drhsqx_dlam = eMatInv * drhsqx_dlam;
            
            qyMatrix = eMatInv * Key;
            %qyrhs = eMatInv * rhsqy;
            drhsqy_dlam = eMatInv * drhsqy_dlam;
            
            %F = rhsu - uqx * qxrhs - uqy * qyrhs;
            dF_dlam = drhsu_dlam - uqx * drhsqx_dlam - uqy * drhsqy_dlam;
            
            dF = uqx * qxMatrix + uqy * qyMatrix + uu;
            
            %u = dF \ F;
            % $$$ [L,U] = lu(dF);
            % $$$ du_dlam = U\(L \ dF_dlam);
            du_dlam = inv(dF) * dF_dlam;
            
            %qx = qxMatrix * u + qxrhs;
            %qy = qyMatrix * u + qyrhs;
            dqx_dlam = qxMatrix * du_dlam + drhsqx_dlam;
            dqy_dlam = qyMatrix * du_dlam + drhsqy_dlam;
        end
        
        
        function A = gen_hdg_matrix(self)
            self.refel = homg.refel (self.Mesh.dim, self.Mesh.order);
            
            [self.SkelInterior2All, self.SkelAll2Interior, self.Bmaps, self.LIFT, self.VtoF] = self.Mesh.generate_skeleton_maps(self.refel);
            
            % variables
            Nfp = self.refel.Nrp ^ (self.refel.dim - 1);
            Nfaces = self.refel.dim * 2;
            Nv = self.refel.Nrp ^ (self.refel.dim);
            K = prod(self.Mesh.nelems);
            
            nx   = self.Mesh.nx;
            ny   = self.Mesh.ny;
            taur = self.Mesh.taur;
            
            maxnnzeros = Nfaces * 2 * Nfp * Nfp * self.Mesh.Ni_faces;
            II = zeros(maxnnzeros,1);
            JJ = zeros(maxnnzeros,1);
            SS = zeros(maxnnzeros,1);
            
            Nfp2 = Nfp * Nfp;
            nnzeros = 0;
            I = eye(Nfp);
            
            for e = 1:K
                % e1 solution
                [du_dlam, dqx_dlam, dqy_dlam] = self.DifflocalSolver(e);
                for f = 1:Nfaces
                    index = (f-1)*Nfp+1:f*Nfp;
                    idxf = self.Mesh.get_skeletal_face_indices(self.refel, e, f);
                    idxf_interior = self.SkelAll2Interior(idxf);
                    
                    % if this is boundary face, ignore
                    if idxf_interior(1) < eps,
                        continue;
                    end
                    
                    Jf = self.Mesh.geometric_factors_face(self.refel,e,f);
                    idxv = self.Mesh.get_discontinuous_face_indices(self.refel, 1, f);
                    
                    
                    % $$$     % construct the residual
                    % $$$     Fhat = Jf .* (refel.Mr * (qx(idxv) * nx(f) + qy(idxv) * ny(f) + ...
                    % $$$                               taur * (u(idxv) - lam)));
                    du_dlam_f  = du_dlam(idxv, index);
                    dqx_dlam_f = dqx_dlam(idxv,index);
                    dqy_dlam_f = dqy_dlam(idxv,index);
                    
                    dFhat_dlam = self.refel.Mr * (diag(Jf .* nx(f)) * dqx_dlam_f +...
                        diag(Jf .* ny(f)) * dqy_dlam_f +...
                        diag(Jf .* taur)  * (du_dlam_f - I));
                    
                    %    res(SkelAll2Interior(idxf)) = res(SkelAll2Interior(idxf)) + Fhat;
                    % self derivative
                    II(nnzeros+1:nnzeros+Nfp2) = repmat(idxf_interior, 1, Nfp);
                    JJ(nnzeros+1:nnzeros+Nfp2) = repmat(idxf_interior', Nfp, 1);
                    SS(nnzeros+1:nnzeros+Nfp2) = dFhat_dlam;
                    
                    nnzeros = nnzeros + Nfp2;
                    
                    for fn = 1:Nfaces
                        indexn = (fn-1)*Nfp+1:fn*Nfp;
                        idxfn = self.Mesh.get_skeletal_face_indices(self.refel, e, fn);
                        idxf_interiorn = self.SkelAll2Interior(idxfn);
                        
                        % derivative wrt f is already counted in the identity I above
                        if (idxf_interiorn(1) < eps) || (fn == f)
                            continue;
                        end
                        
                        du_dlam_fn  = du_dlam(idxv, indexn);
                        dqx_dlam_fn = dqx_dlam(idxv,indexn);
                        dqy_dlam_fn = dqy_dlam(idxv,indexn);
                        
                        dFhat_dlamn = self.refel.Mr * (diag(Jf .* nx(f)) * dqx_dlam_fn +...
                            diag(Jf .* ny(f)) * dqy_dlam_fn +...
                            diag(Jf .* taur)  * du_dlam_fn );
                        
                        % neighbor derivative
                        II(nnzeros+1:nnzeros+Nfp2) = repmat(idxf_interior, 1, Nfp);
                        JJ(nnzeros+1:nnzeros+Nfp2) = repmat(idxf_interiorn', Nfp, 1);
                        SS(nnzeros+1:nnzeros+Nfp2) = dFhat_dlamn;
                        
                        nnzeros = nnzeros + Nfp2;
                    end
                    
                end
            end
            
            A = sparse(II(1:nnzeros),JJ(1:nnzeros),SS(1:nnzeros), self.Mesh.Ni_faces * Nfp, self.Mesh.Ni_faces * Nfp);
            
            self.K = -A;
            self.is_hDG_solve = true;
        end

        function I = skel_to_cg_matrix (self)
            Nfp = self.refel.Nrp ^ (self.refel.dim - 1);
            
            dof = prod(self.Mesh.nelems*self.refel.N + 1);
            hdg_dof = length(self.SkelInterior2All);
            
            num_elem = prod(self.Mesh.nelems);
            
            nx = self.Mesh.nelems(1);
            ny = self.Mesh.nelems(2);
            
            % TODO: need to replace 1.0 with the actual domain length
            hx = 1.0/nx; hy = 1.0/ny;
            fac = hx*hy / (2*(hx+hy));
            
            % add Bdata to skel to get linear dg
            lamAll = zeros(self.Mesh.Ns_faces * Nfp,1);
            lamAll(self.Bmaps) = 0;
            u_skel = zeros(hdg_dof,1);
            I = zeros(dof,hdg_dof);
            u_cg = zeros(dof,1);
            
            for nhdg = 1:hdg_dof
              u_skel(nhdg) = 1;
              lamAll(self.SkelInterior2All) = u_skel;
              u_cg(:) = 0; 
              
              % loop over elements
              for e=1:num_elem
                % TODO: remove hard-coding with 4 here
                for fid=1:4
                    idx = self.Mesh.get_skeletal_face_indices(self.refel, e, fid);
                    [cg_idx, gfid] = self.Mesh.get_continuous_face_indices(self.refel, e, fid);
                    
                    %----Restrict HDG linear skeleton solution to linear CG mesh-------
                    % (Q_0 mu, v)_0 = (mu, I_1 v)_1
                    Jf = self.Mesh.geometric_factors_face(self.refel, e, fid);
                    Md = self.refel.w .* Jf;
                    Mf = self.refel.q1d' * diag(Md) * self.refel.q1d;
                    rhs = Mf * lamAll(idx);  
                    % rhs = Jf .* (self.refel.Mr * lamAll(idx));
                    
                    if self.SkelAll2Interior(idx(1)) < eps
                      if norm(rhs) > 1.e-14
                        error('rhs must be zero here');
                      end
                        u_cg(cg_idx) = u_cg(cg_idx) + fac * rhs;
                    else
                        u_cg(cg_idx) = u_cg(cg_idx) + fac * rhs;
                    end
                    %------------------------------------------------------------------
                end
            end
            
% $$$             bdy_index = ...
% $$$                 self.Mesh.get_boundary_node_indices(self.Mesh.order);
% $$$             u_cg(bdy_index)  = 0; 
            I(:,nhdg) = u_cg;
            u_skel(nhdg) = 0;
            end
            self.I = I;
        end
        
        function [Q] = cg_to_skel_matrix (self)
            Nfp = self.refel.Nrp ^ (self.refel.dim - 1);
            num_elem = prod(self.Mesh.nelems);
            dof = prod(self.Mesh.nelems*self.refel.N + 1);
            hdg_dof = length(self.SkelInterior2All);
            
            lamAll = zeros(self.Mesh.Ns_faces * Nfp,1);
                        
            u_cg = zeros(dof,1);
            Q = zeros(hdg_dof, dof);
            
            for ncg = 1:dof
              u_cg(ncg)  = 1;
              lamAll(:) = 0;
              % loop over elements
              for e=1:num_elem
                % TODO: remove hard-coding with 4 here
                for fid=1:4
                  idx = self.Mesh.get_skeletal_face_indices(self.refel, e, fid);
                  [cg_idx, gfid] = self.Mesh.get_continuous_face_indices(self.refel, e, fid);
                  
                  % Prolongate CG linear solution to HDG linear skeleton
                  lamAll(idx) = u_cg(cg_idx);
                end
              end

              Q(:,ncg) = lamAll(self.SkelInterior2All);
              u_cg(ncg) = 0;
            end
            self.Q = Q;
        end
        
        function [u,qx,qy] = localSolver (self, e1, lam, forcing)
            % disp(e1);
            m     = self.Mesh;
            refel = self.refel;
            
            taur  = m.taur;
            
            LIFT  = self.LIFT;
            VtoF  = self.VtoF;
            
            nx = m.nx;
            ny = m.ny;
            
            % the number of volume point
            Nv = refel.Nrp ^ (refel.dim);
            
            uu = zeros(Nv,Nv);
            uqx = zeros(Nv,Nv);
            uqy = zeros(Nv,Nv);
            
            rhsqx = zeros(Nv,1);
            rhsqy = zeros(Nv,1);
            rhsu  = zeros(Nv,1);
            
            pts = m.element_nodes(e1, refel);
            [Jv, Dv] = m.geometric_factors(refel, pts);
            
            eMat = m.element_mass(e1, refel, Jv);
            eMatInv = inv(eMat);
            
            % number of faces
            Nfaces = refel.dim * 2;
            
            % compute the forcing
            if ( isnumeric(forcing) )
                rhsu = eMat * forcing( ((e1-1)*Nv+1):(e1*Nv) ); % forcing(:,e1);
            else
                rhsu = eMat * forcing(pts);
            end
            % advection stiffness
            [Kex,Key] = m.element_stiffness_advection(e1, refel, Jv, Dv);
            
            uqx = -Kex;
            uqy = -Key;
            % residual for qx and qy equations
            for f = 1:Nfaces %
                idxf = m.get_skeletal_face_indices(refel, e1, f);
                % geometrix factors at gll points on face
                Jf = m.geometric_factors_face(refel,e1,f);
                
                idxv = m.get_discontinuous_face_indices(refel, 1, f);
                
                % residual due to lambda
                rhsfx = -Jf .* (refel.Mr * lam(idxf)) * nx(f);
                rhsfy = -Jf .* (refel.Mr * lam(idxf)) * ny(f);
                % lift to volume residual q equation
                rhsqx(idxv) = rhsqx(idxv) + rhsfx;
                rhsqy(idxv) = rhsqy(idxv) + rhsfy;
                % lift to volume residual u equation
                rhsu(idxv)  = rhsu(idxv) + ...
                    taur * Jf .* (refel.Mr * lam(idxf));
                
                % lift to volume for uu
                bdry =  LIFT(:,:,f) * (diag(Jf) * VtoF(:,:,f));
                bdryx = LIFT(:,:,f) * (diag(Jf) * VtoF(:,:,f)) * nx(f);
                bdryy = LIFT(:,:,f) * (diag(Jf) * VtoF(:,:,f)) * ny(f);
                
                uu  = uu  + taur * bdry;
                uqx = uqx +        bdryx;
                uqy = uqy +        bdryy;
            end
            
            qxMatrix = eMatInv * Kex;
            qxrhs = eMatInv * rhsqx;
            
            qyMatrix = eMatInv * Key;
            qyrhs = eMatInv * rhsqy;
            
            F = rhsu - uqx * qxrhs - uqy * qyrhs;
            dF = uqx * qxMatrix + uqy * qyMatrix + uu;
            
            u = dF \ F;
            qx = qxMatrix * u + qxrhs;
            qy = qyMatrix * u + qyrhs;
            
        end
        
        function u_cg = plot_hdg_skel (self, u_hat)
          % ignores BData
          
          Nfp = self.refel.Nrp ^ (self.refel.dim - 1);
          num_elem = prod(self.Mesh.nelems);
          
          lamAll = zeros(self.Mesh.Ns_faces * Nfp,1);
          dof = prod(self.Mesh.nelems*self.refel.N + 1);
          
          % create zero array of CG size
          u_cg = zeros(dof,1);

          lamAll(self.SkelInterior2All) = u_hat; 

          % loop over elements
          for e=1:num_elem
            % TODO: remove hard-coding with 4 here
            for fid=1:4
              idx = self.Mesh.get_skeletal_face_indices(self.refel, e, fid);
              [cg_idx, gfid] = self.Mesh.get_continuous_face_indices(self.refel, e, fid);
              
              % Prolongate CG linear solution to HDG linear skeleton
              u_cg(cg_idx) = u_cg(cg_idx) + lamAll(idx);
            end
          end
          
          u_cg = reshape(u_cg, self.Mesh.nelems(1)+1, self.Mesh.nelems(2)+1);
          % scale values ...
          u_cg(2:end-1,2:end-1) = u_cg(2:end-1,2:end-1) * 0.125;
          
          u_cg(1,2:end-1) = u_cg(1,2:end-1)*0.5;
          u_cg(end,2:end-1) = u_cg(end,2:end-1)*0.5;
          
          u_cg(2:end-1,1) = u_cg(2:end-1,1)*0.5;
          u_cg(2:end-1,end) = u_cg(2:end-1,end)*0.5;
          
          % draw
          subplot(1,2,1);
          imagesc(u_cg);
          colorbar;

          % special plot ...
          gamma = zeros(self.Mesh.nelems(1)*Nfp + self.Mesh.nelems(1) + 1, ...
                        self.Mesh.nelems(2)*Nfp + self.Mesh.nelems(2) + 1);
          
          for e=1:num_elem
            for fid=1:4
              idx = self.Mesh.get_skeletal_face_indices(self.refel, e, fid);
              [i,j] = ind2sub (self.Mesh.nelems, e);
              
              switch fid
                case 1
                  xidx = (i-1)*(Nfp+1) + 1;
                  yidx = (j-1)*(Nfp+1) + 1 + (1:Nfp);
                  
                case 2
                  xidx = i*(Nfp+1) + 1;
                  yidx = (j-1)*(Nfp+1) + 1 + (1:Nfp);
                  
                case 3
                  xidx = (i-1)*(Nfp+1) + 1 + (1:Nfp);
                  yidx = (j-1)*(Nfp+1) + 1;
                case 4
                  xidx = (i-1)*(Nfp+1) + 1 + (1:Nfp);
                  yidx = j*(Nfp+1) + 1;
              end
              
              gamma(xidx, yidx) = lamAll(idx);
            end
          end
          subplot(1,2,2);
          imagesc(gamma); colorbar;
        end % plot_hdg_skel
        
        function u_hat_z = clear_skel_boundary(self, u_hat)
          Nfp = self.refel.Nrp ^ (self.refel.dim - 1);
          
          Ni = self.Mesh.Ni_faces;
          nelems = self.Mesh.nelems;
          
          nxf = nelems(2)*(nelems(1)+1);
          
          lamAll = zeros(self.Mesh.Ns_faces * Nfp,1);
          
          lamAll(self.SkelInterior2All) = u_hat;
          
          % clear lower x-faces
          for i=1:(nelems(1)+1)
            lamAll((i-1)*Nfp + 1) = 0;
          end
          %clear upper x-faces
          for i=((nelems(1)+1)*(nelems(2)-1)+1):(nelems(2)*(nelems(1)+1))
            lamAll(i*Nfp) = 0;
          end
          
%           % clear lower x-faces
          for i=1:(nelems(2)+1)
            lamAll((nxf+nelems(1))*Nfp + (i-1)*nelems(1)*Nfp + 1) = 0;
            lamAll((nxf+nelems(1))*Nfp + (i-1)*nelems(1)*Nfp) = 0;
          end
%           %clear upper x-faces
%           for i=((nelems(2)+1)*(nelems(1)-1)+1):(nelems(1)*(nelems(2)+1))
%             lamAll((nxf)*Nfp + (i-1)*nelems(1)*Nfp) = 6;
%           end

          u_hat_z = lamAll(self.SkelInterior2All);
          
        end
        
    end %methods
    
end %classdef

