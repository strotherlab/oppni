function out = ANOVA2( datamat, design )
%
% . 2-way ANOVA --> 2 experimental manipulation, multiple levels
%

% signif = 0.05;


if( size(design,2) <2 )
    error('you only have 1 design factor -- run 1-way ANOVA instead!');
elseif( size(design,2)  >2 )
    error('can only have 2 design factors');
end


NF=2;
for(f=1:2) NL(f) = numel( unique(design(:,f)) ); end
for(i=1:NL(1))
for(j=1:NL(2))
    rr(i,j) = sum( design(:,1)==i & design(:,2)==j );
end
end
for(f=1:2) 
for(i=1:NL(f))
    rf{f}(i)= sum( design(:,f)==i );
end
end

ssw = 0;
for(i=1:NL(1))
for(j=1:NL(2))

	yijk = datamat(:, design(:,1)==i & design(:,2)==j );
	yij  = mean(yijk,2);
	ssw  = ssw+ sum( (yijk-mean(yij)).^2,2 );
end
end
dfw = sum(rr(:)) - prod(NL);

y_gm = mean(datamat,2);

ssg1=0;
for(i=1:NL(1))
	ssg1 = ssg1 + rf{1}(i)*( mean(datamat(:,design(:,1)==i),2) - mean(datamat,2) ).^2;
end
dfg1 = NL(1)-1;

ssg2=0;
for(j=1:NL(2))
	ssg2 = ssg2 + rf{2}(j)*( mean(datamat(:,design(:,2)==j),2) - mean(datamat,2) ).^2;
end
dfg2 = NL(2)-1;

ssi=0;
for(i=1:NL(1))
for(j=1:NL(2))
	yij = mean( datamat( :, design(:,1)==i & design(:,2)==j ),2 );
	yi  = mean( datamat( :, design(:,1)==i ),2 );
	yj  = mean( datamat( :, design(:,2)==j ),2 );

	ssi = ssi + rr(i,j)*( yij - yi - yj +y_gm ).^2;
end
end
dfi = (NL(1)-1)*(NL(2)-1);

ms=[ ssg1./dfg1, ssg2./dfg2, ssi./dfi, ssw./dfw ];
ff=bsxfun(@rdivide, ms(:,1:3), ms(:,4) );
PP = [ 1-fcdf(ff(:,1),dfg1,dfw), 1-fcdf(ff(:,2),dfg2,dfw), 1-fcdf(ff(:,3),dfi,dfw)];

%% statistics - treat1 / treat2 / interaction

% group1 effect
out.fstat_G1    = ff(:,1);
out.fstat_G1_p  = PP(:,1);
% group2 effect
out.fstat_G2    = ff(:,2);
out.fstat_G2_p  = PP(:,2);
% interaction effect
out.fstat_G12   = ff(:,3);
out.fstat_G12_p = PP(:,3);


out.testname = 'anova2';
