classdef xform
    %XFORM transformation functions
    %   of the form, Xout = foo(Xin)
    
    methods(Static)
        function Xout = identity (Xin)
            Xout = Xin;
        end
        
        function Xout = twoX (Xin)
            Xout = Xin;
            Xout(:,1) = 2*Xout(:,1);
        end
        
        function Xout = twoY (Xin)
            Xout = Xin;
            Xout(:,2) = 2*Xout(:,2);
        end
        
        function Xout = twoSix (Xin)
            % currently only 2D
            Xout = 2*Xin;
            Xout(:,2) = 3*Xout(:,2);
        end
        
        function Xout = bowl (Xin)
            Xtemp(:,1) = 2*Xin(:,1) - 1;
            Xtemp(:,2) = 2*Xin(:,2) - 1;
            Xtemp(:,3) = Xin(:,3) - 1;
            tempX = sqrt(1.0 - 0.5*Xtemp(:,2).^2 - 0.5*Xtemp(:,3).^2 + 0.3333*Xtemp(:,2).^2.*Xtemp(:,3).^2);
            tempY = sqrt(1.0 - 0.5*Xtemp(:,3).^2 - 0.5*Xtemp(:,1).^2 + 0.3333*Xtemp(:,3).^2.*Xtemp(:,1).^2);
            tempZ = sqrt(1.0 - 0.5*Xtemp(:,1).^2 - 0.5*Xtemp(:,2).^2 + 0.3333*Xtemp(:,1).^2.*Xtemp(:,2).^2);
            Xout(:,1) = Xtemp(:,1).*tempX;
            Xout(:,2) = Xtemp(:,2).*tempY;
            Xout(:,3) = Xtemp(:,3).*tempZ;
        end
        
        function Xout = pcylhalf (Xin)
            Xtemp(:,1) = 2*Xin(:,1) - 1;
            Xtemp(:,2) = 2*Xin(:,2) - 1;
            Xtemp(:,3) = Xin(:,3) - 1;
            R = -125/27*Xtemp(:,3).^3 - 125/18*Xtemp(:,3).^2 - 20/9*Xtemp(:,3) + 43/54;
            I = find(-1*Xtemp(:,3) < 0.2);
            R(I) = ones(size(I));
            I = find(-1*Xtemp(:,3) > 0.8);
            R(I) = 0.5*ones(size(I));
            tempX = R.*sqrt(1 - 0.5*Xtemp(:,2).^2);
            tempY = R.*sqrt(1 - 0.5*Xtemp(:,1).^2);
            Xout(:,1) = Xtemp(:,1).*tempX;
            Xout(:,2) = Xtemp(:,2).*tempY;
            Xout(:,3) = Xtemp(:,3);
        end
        
        function Xout = pcyldouble (Xin)
            Xtemp(:,1) = 2*Xin(:,1) - 1;
            Xtemp(:,2) = 2*Xin(:,2) - 1;
            Xtemp(:,3) = Xin(:,3) - 1;
            R = 250/27*Xtemp(:,3).^3 + 125/9*Xtemp(:,3).^2 + 40/9*Xtemp(:,3) + 38/27;
            I = find(-1*Xtemp(:,3) < 0.2);
            R(I) = ones(size(I));
            I = find(-1*Xtemp(:,3) > 0.8);
            R(I) = 2.0*ones(size(I));
            tempX = R.*sqrt(1 - 0.5*Xtemp(:,2).^2);
            tempY = R.*sqrt(1 - 0.5*Xtemp(:,1).^2);
            Xout(:,1) = Xtemp(:,1).*tempX;
            Xout(:,2) = Xtemp(:,2).*tempY;
            Xout(:,3) = Xtemp(:,3);
        end
        
        function Xout = disc3 (Xin)
            d = size(Xin,2);
            if (d == 2)
                Xtemp(:,1) = 2*Xin(:,1) - 1;
                Xtemp(:,2) = 2*Xin(:,2) - 1;
                I = find(Xtemp(:,1).^2+Xtemp(:,2).^2);
                notI = find(not(Xtemp(:,1).^2+Xtemp(:,2).^2));
                Xtemp1 = (zeros([length(Xtemp),1]));
                Xtemp1(I) = sqrt(Xtemp(I,1).^2 + Xtemp(I,2).^2 - Xtemp(I,1).^2.*Xtemp(I,2).^2);
                Xtemp1(I) = Xtemp1(I) ./ sqrt(Xtemp(I,1).^2 + Xtemp(I,2).^2);
                Xtemp1(notI) = 0;
                Xtemp1(notI) = 0;
                Xout = zeros([length(Xtemp),2]);
                Xout(:,1) = Xtemp(:,1) .* Xtemp1;
                Xout(:,2) = Xtemp(:,2) .* Xtemp1;
            end
        end
        
        function Xout = cylinder (Xin)
            Xtemp(:,1) = 2*Xin(:,1) - 1;
            Xtemp(:,2) = 2*Xin(:,2) - 1;
            I = find(Xtemp(:,1).^2+Xtemp(:,2).^2);
            notI = find(not(Xtemp(:,1).^2+Xtemp(:,2).^2));
            Xtemp1 = (zeros([length(Xtemp),1]));
            Xtemp1(I) = sqrt(Xtemp(I,1).^2 + Xtemp(I,2).^2 - Xtemp(I,1).^2.*Xtemp(I,2).^2);
            Xtemp1(I) = Xtemp1(I) ./ sqrt(Xtemp(I,1).^2 + Xtemp(I,2).^2);
            Xtemp1(notI) = 0;
            Xtemp1(notI) = 0;
            Xout = zeros([length(Xtemp),2]);
            Xout(:,1) = Xtemp(:,1) .* Xtemp1;
            Xout(:,2) = Xtemp(:,2) .* Xtemp1;
            Xout(:,3) = Xin(:,3);
        end
        
        function Xout = shell (Xin)
            d = size(Xin, 2);
            R2 = 1.0; % hard coded for now.
            R1 = 0.7; % hard coded for now.
            R2byR1 = R2 / R1;
            R1sqrbyR2 = R1 * R1 / R2;
            
            if (d == 2)
                x = zeros( size(Xin(:,1)) );
                y = tan ( Xin(:,1)  * pi/4 );
                R = R1sqrbyR2 * ( R2byR1.^(Xin(:,2) + 1) ) ;
            else
                x = tan ( Xin(:,1)  * pi/4 );
                y = tan ( Xin(:,2)  * pi/4 );
                R = R1sqrbyR2 * ( R2byR1.^(Xin(:,3) + 1) );
            end
            
            q = R ./ sqrt (x.*x + y.*y + 1);
            
            if (d == 3)
                Xout(:,1) =  q.* y;
                Xout(:,2) = -q.* x;
                Xout(:,3) =  q;
            else
                Xout(:,1) =   q.* y;
                Xout(:,2) =   q;
            end
        end
    end % Static methods
    
end

% for (il = 0; il < np; ++il) {
%       /* transform abc[0] and y in-place for nicer grading */
%       x = tan (EX[il] * M_PI_4);
%       y = tan (EY[il] * M_PI_4);
%
%       /* compute transformation ingredients */
%       R = R1sqrbyR2 * pow (R2byR1, EZ[il]);
%       q = R / sqrt (x * x + y * y + 1.);
%
%       /* assign correct coordinates based on patch id */
%       X[il] = +q * y;
%       Y[il] = -q * x;
%       Z[il] = +q;
%     }
