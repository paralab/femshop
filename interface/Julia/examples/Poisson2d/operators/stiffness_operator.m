function [K] = stiffness_operator(u, args)
    if args.left==true
        K=args.mesh.assemble_stiffness(args.config.basis_order_min);
    else
        K=args.mesh.assemble_stiffness(args.config.basis_order_min) * u;
    end
end