function [M] = mass_operator(u, args)
    if args.left==true
        M=args.mesh.assemble_mass(args.config.basis_order_min);
    else
        M=args.mesh.assemble_mass(args.config.basis_order_min) * u;
    end
end