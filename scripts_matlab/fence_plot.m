function fence_plot( colData, refchange, tolerance )
%
% This script gives a "fence plot", as alternative to
% boxplot --> individual sample changes are plotted
%
%        fence_plot( colData, refchange, tolerance )
%

newfig=0;
[Nsamp Ntreat] = size(colData);

%
if(newfig>0) figure; hold on;
else                 hold on;
end

if(newfig>1)  
    boxplot(colData,'whisker',0,'symbol','.k','color','k');
end
    

for( n=1:Ntreat )
    
% %     plot( n*ones(Nsamp,1), colData(:,n ), 'ok' );
    
    if( refchange == 0 )    
        if(n<Ntreat)
        for(m=1:Nsamp)
           if( colData(m,n+1) > (colData(m,n)    + tolerance) )
              plot( [n n+1], [colData(m,n) colData(m,n+1)], '.-b', 'linewidth',2 );
              plot( [n n+1], [colData(m,n) colData(m,n+1)], 'ok', 'linewidth',1, 'markersize',4 );
           elseif( colData(m,n+1) < (colData(m,n) - tolerance) )
              plot( [n n+1], [colData(m,n) colData(m,n+1)], '.-r', 'linewidth',2 );            
              plot( [n n+1], [colData(m,n) colData(m,n+1)], 'ok', 'linewidth',1, 'markersize',4 );
           else
              plot( [n n+1], [colData(m,n) colData(m,n+1)], ':k', 'linewidth',1 );
           end
        end
        end
    else

        for(m=1:Nsamp)
           if(n<Ntreat)           
            plot( [n n+1], [colData(m,n) colData(m,n+1)], '-k' );
           end
           if( colData(m,n) > ( colData(m,refchange) + tolerance ) )
                plot( n, colData(m,n), 'ob' );
           elseif(colData(m,n) < ( colData(m,refchange) - tolerance ) )
                plot( n, colData(m,n), 'or' );
           end
        end
    end
    
end

xlim([0.5 Ntreat+0.5]);
