#=
# Applies various operators
=#

export dg_mass_inv_advective
export dg_surface_int

# M_inv * Ad * u
function dg_mass_inv_advective(u)
    return rx.*(refel.Dr*u);
end

# surface integral(normal . u * v)
function dg_surface_int(u)
    return refel.lift*(Fscale.*(normals.*u));
end