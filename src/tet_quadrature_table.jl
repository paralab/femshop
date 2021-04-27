#=
Table of quadrature points and weights for the reference tet.

Formulae from:

Hermann Engels,
Numerical Quadrature and Cubature,
Academic Press, 1980,
ISBN: 012238850X,
LC: QA299.3E5. 

=#

function tetrahedron_quadrature_nodes_weights(order)
    
    if order > 15 || order < 1
        printerr("unsupported order for tet quadrature nodes. Choose 1 to 15.")
        return zeros(0,0);
    end
    
    cubTet = Array{Array{Float64,2},1}(undef,15);
    
    cubTet[1] = [-0.5 -0.5 -0.5 2.0];
    
    cubTet[2] = [
        0.17082039324993703 -0.7236067977499789 -0.7236067977499789 0.25; 
        -0.7236067977499789 -0.7236067977499789 -0.7236067977499789 0.25; 
        -0.7236067977499789 -0.7236067977499789 0.17082039324993703 0.25; 
        -0.7236067977499789 0.17082039324993703 -0.7236067977499789 0.25
        ];
    
    cubTet[3] = cubTet[2];
    cubTet[4] = cubTet[2];
    
    cubTet[5] = [
        -0.5 -0.5 -0.5 -0.8; 
        0.0 -0.6666666666666665 -0.6666666666666665 0.45; 
        -0.6666666666666665 -0.6666666666666665 -0.6666666666666665 0.45; 
        -0.6666666666666665 -0.6666666666666665 0.0 0.45; 
        -0.6666666666666665 0.0 -0.6666666666666665 0.45
        ];
    
    cubTet[6] = [
        0.1368611683936889 -0.7122870561312296 -0.7122870561312296 0.2177650698804054; 
        -0.7122870561312296 -0.7122870561312296 -0.7122870561312296 0.2177650698804054; 
        -0.7122870561312296 -0.7122870561312296 0.1368611683936889 0.2177650698804054; 
        -0.7122870561312296 0.1368611683936889 -0.7122870561312296 0.2177650698804054; 
        -1.0 0.0 0.0 0.0214899534130631; 
        0.0 -1.0 0.0 0.0214899534130631; 
        0.0 0.0 -1.0 0.0214899534130631; 
        0.0 -1.0 -1.0 0.0214899534130631; 
        -1.0 0.0 -1.0 0.0214899534130631; 
        -1.0 -1.0 0.0 0.0214899534130631
        ];
    
    cubTet[7] = cubTet[6];
    cubTet[8] = cubTet[6];
    cubTet[9] = cubTet[6];
    cubTet[10] = cubTet[6];
    
    cubTet[11] = [
        -0.5 -0.5 -0.5 0.1817020685825351; 
        -1.0 -0.33333333333333337 -0.33333333333333337 0.0361607142857143; 
        -0.33333333333333337 -0.33333333333333337 -0.33333333333333337 0.0361607142857143; 
        -0.33333333333333337 -0.33333333333333337 -1.0 0.0361607142857143; 
        -0.33333333333333337 -1.0 -0.33333333333333337 0.0361607142857143; 
        0.4545454545454546 -0.8181818181818182 -0.8181818181818182 0.0698714945161738; 
        -0.8181818181818182 -0.8181818181818182 -0.8181818181818182 0.0698714945161738; 
        -0.8181818181818182 -0.8181818181818182 0.4545454545454546 0.0698714945161738; 
        -0.8181818181818182 0.4545454545454546 -0.8181818181818182 0.0698714945161738; 
        -0.13310030714732857 -0.8668996928526714 -0.8668996928526714 0.0656948493683187; 
        -0.8668996928526714 -0.13310030714732857 -0.8668996928526714 0.0656948493683187; 
        -0.8668996928526714 -0.8668996928526714 -0.13310030714732857 0.0656948493683187; 
        -0.8668996928526714 -0.13310030714732857 -0.13310030714732857 0.0656948493683187; 
        -0.13310030714732857 -0.8668996928526714 -0.13310030714732857 0.0656948493683187; 
        -0.13310030714732857 -0.13310030714732857 -0.8668996928526714 0.0656948493683187
        ];
    
    cubTet[12] = cubTet[11];
    cubTet[13] = cubTet[11];
    cubTet[14] = cubTet[11];
    cubTet[15] = cubTet[11];
    
    return cubTet[order];
end