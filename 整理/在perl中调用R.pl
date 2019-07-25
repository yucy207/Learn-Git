use Statistics::R;

    my $coveragelist=join ",",@covers;
    my $samplenamelist=join ",",@samples;
    my $samplenum=@samples;
    my $pic="$imgdir/mapping_coverage.png";
    my $R = Statistics::R->new();
    $R->startR;
    $R->send(qq` myColor=c("#aa4643","#89a54e","#71588f","#4198af","#db843d") `);
    $R->send(qq` png(\"$pic\",width=800,height=600) `) ;
    $R->send(qq` barplot(c($coveragelist),col=myColor[1],names.arg=c($samplenamelist),xlab="Samples",ylab="Percentage",lwd=2) `) if($samplenum<=5);
    $R->send(qq` barplot(c($coveragelist),col=myColor[1],xlab="Samples",ylab="Percentage",lty=1,lwd=2) `) if($samplenum>5);
    $R->send(qq` legend("topleft",c("MEAN_TARGET_COVERAGE"),cex=1,lwd=2,col=myColor[1],lty=c(1)) `);
    $R->send(qq` dev.off() `);  
}