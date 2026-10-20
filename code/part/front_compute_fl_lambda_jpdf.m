function result = front_compute_fl_lambda_jpdf(fronts, data, binEdges)
%FRONT_COMPUTE_FL_LAMBDA_JPDF Diagnose every retained FL point on lambda-JPDF.

if nargin < 3 || isempty(binEdges), binEdges = -5:0.05:5; end
required = {'lonir','latir','pupx','pupy','pvpx','pvpy','pbhpx','pbhpy','f'};
for k=1:numel(required)
    if ~isfield(data,required{k})
        error('Front:MissingVariable','Required field missing: %s',required{k});
    end
end

flId=[]; pointId=[]; x=[]; y=[];
for k=1:numel(fronts)
    n=numel(fronts(k).ridge.x);
    flId=[flId;repmat(fronts(k).fl_id,n,1)]; %#ok<AGROW>
    pointId=[pointId;(1:n).']; %#ok<AGROW>
    x=[x;fronts(k).ridge.x(:)]; %#ok<AGROW>
    y=[y;fronts(k).ridge.y(:)]; %#ok<AGROW>
end

ux=interp2(data.lonir,data.latir,data.pupx,x,y,'linear');
uy=interp2(data.lonir,data.latir,data.pupy,x,y,'linear');
vx=interp2(data.lonir,data.latir,data.pvpx,x,y,'linear');
vy=interp2(data.lonir,data.latir,data.pvpy,x,y,'linear');
bx=interp2(data.lonir,data.latir,data.pbhpx,x,y,'linear');
by=interp2(data.lonir,data.latir,data.pbhpy,x,y,'linear');
ff=interp2(data.lonir,data.latir,data.f,x,y,'linear');
n=numel(x);
lambdaMinus=nan(n,1); lambdaPlus=nan(n,1);
theta=nan(n,1); epsilon=nan(n,1); className=strings(n,1);
deltaF=nan(n,1); lambdaJpdfY=nan(n,1); valid=false(n,1);
tol=1e-10;
for k=1:n
    z=[ux(k) uy(k) vx(k) vy(k) bx(k) by(k) ff(k)];
    if any(~isfinite(z)) || ff(k)==0 || hypot(bx(k),by(k))<=eps, continue; end
    A=[ux(k) uy(k);vx(k) vy(k)];
    [V,D]=eig(A.'); ev=diag(D);
    if max(abs(imag(ev)))>tol, continue; end
    [ev,ord]=sort(real(ev),'ascend'); V=real(V(:,ord));
    if abs(ev(2)-ev(1))<=100*eps(max(1,norm(A,'fro'))), continue; end
    sig=hypot(ux(k)-vy(k),uy(k)+vx(k)); zet=vx(k)-uy(k);
    omega=sig.^2-zet.^2;
    if omega < -tol*max(1,sig.^2+zet.^2) || sig<=eps, continue; end
    omega=max(omega,0);
    lambdaMinus(k)=ev(1); lambdaPlus(k)=ev(2);
    deltaF(k)=(ux(k)+vy(k))/ff(k);
    lambdaJpdfY(k)=-sqrt(omega)/abs(ff(k));
    theta(k)=acosd(max(-1,min(1,-zet/sig)));
    vminus=V(:,1)/norm(V(:,1)); gradb=[bx(k);by(k)]; gradb=gradb/norm(gradb);
    epsilon(k)=acosd(max(0,min(1,abs(dot(vminus,gradb)))));
    if ev(1)>0
        className(k)="Ds+";
    elseif ev(2)<0
        className(k)="Cs-";
    elseif ev(2)>0 && abs(ev(2))>abs(ev(1))
        className(k)="Sd-";
    elseif ev(2)>0 && abs(ev(1))>abs(ev(2))
        className(k)="Sc-";
    else
        className(k)="boundary";
    end
    valid(k)=true;
end

realValid=valid;
inRange=realValid & deltaF>=binEdges(1) & deltaF<=binEdges(end) & ...
    lambdaJpdfY>=binEdges(1) & lambdaJpdfY<=binEdges(end);
[counts,~,~,ix,iy]=histcounts2(deltaF(inRange),lambdaJpdfY(inRange),binEdges,binEdges);
thetaValid=theta(inRange); epsilonValid=epsilon(inRange); inside=ix>0 & iy>0;
sz=[numel(binEdges)-1 numel(binEdges)-1]; subs=[ix(inside) iy(inside)];
thetaMean=accumarray(subs,thetaValid(inside),sz,@mean,NaN);
thetaMedian=accumarray(subs,thetaValid(inside),sz,@median,NaN);
epsilonMean=accumarray(subs,epsilonValid(inside),sz,@mean,NaN);
epsilonMedian=accumarray(subs,epsilonValid(inside),sz,@median,NaN);
bw=diff(binEdges(1:2)); pdf=counts/(sum(counts(:))*bw^2); pdf(counts==0)=NaN;

points=table(flId,pointId,x,y,deltaF,lambdaJpdfY,lambdaMinus,lambdaPlus, ...
    className,theta,epsilon,realValid,inRange,'VariableNames',{'fl_id','fl_point_index', ...
    'x_km','y_km','delta_over_f','lambda_jpdf_y','AT_lambda_minus_s', ...
    'AT_lambda_plus_s','lambda_class','theta_deg','epsilon_deg','real_eigen_valid','in_jpdf_range'});
centers=0.5*(binEdges(1:end-1)+binEdges(2:end));
result=struct('points',points,'bin_edges',binEdges,'bin_centers',centers, ...
    'counts',counts,'pdf',pdf,'theta_mean',thetaMean,'theta_median',thetaMedian, ...
    'epsilon_mean',epsilonMean,'epsilon_median',epsilonMedian,'sample_count',nnz(realValid), ...
    'in_range_count',nnz(inRange));
end

