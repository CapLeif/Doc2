function [lambdaMinus, lambdaPlus, valid] = front_compute_at_lambda_minus(x, y, data)
%FRONT_COMPUTE_AT_LAMBDA_MINUS Evaluate the real hyperbolic A' eigenpair.

ux = interp2(data.lonir,data.latir,data.pupx,x,y,'linear');
uy = interp2(data.lonir,data.latir,data.pupy,x,y,'linear');
vx = interp2(data.lonir,data.latir,data.pvpx,x,y,'linear');
vy = interp2(data.lonir,data.latir,data.pvpy,x,y,'linear');

discriminant = (ux-vy).^2 + 4.*uy.*vx;
tolerance = 100 .* eps(max(1,abs(ux)+abs(uy)+abs(vx)+abs(vy))).^2;
realPair = isfinite(discriminant) & discriminant >= -tolerance;
rootValue = sqrt(max(discriminant,0));
lambdaMinus = 0.5 .* (ux+vy-rootValue);
lambdaPlus = 0.5 .* (ux+vy+rootValue);
separated = rootValue > 100 .* eps(max(1,abs(ux)+abs(uy)+abs(vx)+abs(vy)));
valid = realPair & separated & isfinite(lambdaMinus) & isfinite(lambdaPlus) & ...
    lambdaMinus < 0;
lambdaMinus(~valid) = NaN;
lambdaPlus(~valid) = NaN;
end

