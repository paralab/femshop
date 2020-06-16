#=
# Applies various operators
=#

export mass_inv_advective
export surface_int

# M_inv * Ad * u
function mass_inv_advective(u)
    return rx.*(refel.Dr*u);
end

# surface integral(normal . u * v)
function surface_int(u)
    return refel.lift*(Fscale.*(normals.*u));
end