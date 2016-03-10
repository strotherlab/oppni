function out = module_ANOVA2( datamat, des, signif )

Ymat = datamat;

NF=2;
for(f=1:2) NL(f) = numel( unique(des(:,f)) ); end
for(i=1:NL(1))
for(j=1:NL(2))
    rr(i,j) = sum( des(:,1)==i & des(:,2)==j );
end
end
for(f=1:2) 
for(i=1:NL(f))
    rf{f}(i)= sum( des(:,f)==i );
end
end

ssw = 0;
for(i=1:NL(1))
for(j=1:NL(2))

	yijk = Ymat(:, des(:,1)==i & des(:,2)==j );
	yij  = mean(yijk,2);
	ssw  = ssw+ sum( (yijk-mean(yij)).^2,2 );
end
end
dfw = sum(rr(:)) - prod(NL);

y_gm = mean(Ymat,2);

ssg1=0;
for(i=1:NL(1))
	ssg1 = ssg1 + rf{1}(i)*( mean(Ymat(:,des(:,1)==i),2) - mean(Ymat,2) ).^2;
end
dfg1 = NL(1)-1;

ssg2=0;
for(j=1:NL(2))
	ssg2 = ssg2 + rf{2}(j)*( mean(Ymat(:,des(:,2)==j),2) - mean(Ymat,2) ).^2;
end
dfg2 = NL(2)-1;

ssi=0;
for(i=1:NL(1))
for(j=1:NL(2))
	yij = mean( Ymat( :, des(:,1)==i & des(:,2)==j ),2 );
	yi  = mean( Ymat( :, des(:,1)==i ),2 );
	yj  = mean( Ymat( :, des(:,2)==j ),2 );

	ssi = ssi + rr(i,j)*( yij - yi - yj +y_gm ).^2;
end
end
dfi = (NL(1)-1)*(NL(2)-1);

ms=[ ssg1./dfg1, ssg2./dfg2, ssi./dfi, ssw./dfw ];
ff=bsxfun(@rdivide, ms(:,1:3), ms(:,4) );
PP = [ 1-fcdf(ff(:,1),dfg1,dfw), 1-fcdf(ff(:,2),dfg2,dfw), 1-fcdf(ff(:,3),dfi,dfw)];

out.p_f = PP;
out.Fstat = ff;