package package::gene;
use strict;
use warnings;

sub run{
    my $GENE       = shift @_;
    my $report_dir = shift @_;
    my $name       = shift @_;
    my $format     = shift @_;
    my $library    = shift @_;
    my $annotation = shift @_;
    my $analysis   = shift @_;
    my $filter     = shift @_;
    my $run        = shift @_;
    my $DoSigGene  = (exists $run->{"DoSigGene"}{1}) ? $run->{"DoSigGene"}{1} : "FALSE";
    my $significant_gene = ();
    if($DoSigGene eq "TRUE"){
        mkdir $report_dir."/document/2_SNV/significant_gene" if(not -e $report_dir."/document/2_SNV/significant_gene");
        $significant_gene = $report_dir."/document/2_SNV/significant_gene/$name.sig_gene.txt";
    }
    $GENE->set_row(0, 60);
    my @case=();@case=split /,/,$run->{"Case"}{1} if(exists $run->{"Case"}{1});
    my @control=();@control=split /,/,$run->{"Control"}{1} if(exists $run->{"Control"}{1});
    my @sample=(@case,@control);
    my @array = ("GeneskyExonDB SNV Count","GeneskyExonDB Mutation(0|1|2)","OMIM","HPO","HGMD_gene","MGI","GO_BP","GO_MF","GO_CC","KEGG_Pathway","Gene_Info");
    my @title=("Gene","First Priority","SNV count","Sample count","Mutation(0|1|2)","p1","p2","p3","p4","Mutation 1","Mutation 2",@array);
    my @title2=("Gene","First Priority","SNV count","Sample count","Mutation(0|1|2)","Case(0|1|2)","Control(0|1|2)","p1","p2","p3","p4","Mutation 1","Mutation 2",@array);
    my $row=0;
    if(exists($control[0])){
        for my $i(0..@title2-1){
            $GENE->write($row,$i,$title2[$i],$format->{"title"});
        }
    }else{
        for my $i(0..@title-1){
            $GENE->write($row,$i,$title[$i],$format->{"title"});
        }
    }

    $row++;
    my %data=();
    my %casehash=();
    my %controlhash=();
    foreach (@case){
        $casehash{$_}++;
    }
    foreach (@control){
        $controlhash{$_}++;
    }
    
    my @control_model = ("normal","strict","normal_pass","strict_pass");
    foreach my $id(sort{$a <=> $b}keys %{$filter})
    {
        foreach my $title(keys %{$filter->{$id}})
        {
            foreach my $alt(keys %{$filter->{$id}{$title}})
            {
                my $Gene=$library->{$title}{$alt}{"Gene"};
                $data{$Gene}{"Gene"}=$Gene;
                my %mut=();
                foreach(split /,/,$analysis->{$title}{$alt}{"Mutation 0"}){$mut{(split /:/,$_)[0]}=0 if($_=~ /\w/);}
                foreach(split /,/,$analysis->{$title}{$alt}{"Mutation 1"}){$mut{(split /:/,$_)[0]}=1 if($_=~ /\w/);}
                foreach(split /,/,$analysis->{$title}{$alt}{"Mutation 2"}){$mut{(split /:/,$_)[0]}=2 if($_=~ /\w/);}
                for my $i(0..@sample-1){
                    next if(not exists $mut{$sample[$i]});
                    $data{$Gene}{"Sample"}{$sample[$i]}+=$mut{$sample[$i]};
                    $data{$Gene}{"Sample"}{$sample[$i]}=2 if($data{$Gene}{"Sample"}{$sample[$i]}>2);
                    $data{$Gene}{"Site"}{$title}=1;
                    $data{$Gene}{"Priority"}{$sample[$i]}{$filter->{$id}{$title}{$alt}{"First Priority"}}+=$mut{$sample[$i]};
                    foreach(@array){
                        next if($_ eq "GeneskyExonDB SNV Count" or $_ eq "GeneskyExonDB Mutation(0|1|2)");
                        $data{$Gene}{$_} = $annotation->{$title}{$alt}{$_} if(exists $annotation->{$title}{$alt}{$_});
                    }
                    $data{$Gene}{"GeneskyExonDB SNV Count"}       = $annotation->{$title}{$alt}{"GeneskyExonDB SNV Count(normal)"} if(exists $annotation->{$title}{$alt}{"GeneskyExonDB SNV Count(normal)"});
                    $data{$Gene}{"GeneskyExonDB Mutation(0|1|2)"} = $annotation->{$title}{$alt}{"GeneskyExonDB Mutation(0|1|2)(normal)"} if(exists $annotation->{$title}{$alt}{"GeneskyExonDB Mutation(0|1|2)(normal)"});
                    foreach(@control_model){
                        $data{$Gene}{"GeneskyExonDB SNV Count"}       = $annotation->{$title}{$alt}{"GeneskyExonDB SNV Count($_)"} if(exists $run->{"GeneControlModel"}{1} and $run->{"GeneControlModel"}{1} eq $_ and exists $annotation->{$title}{$alt}{"GeneskyExonDB SNV Count($_)"});
                        $data{$Gene}{"GeneskyExonDB Mutation(0|1|2)"} = $annotation->{$title}{$alt}{"GeneskyExonDB Mutation(0|1|2)($_)"} if(exists $run->{"GeneControlModel"}{1} and $run->{"GeneControlModel"}{1} eq $_ and exists $annotation->{$title}{$alt}{"GeneskyExonDB Mutation(0|1|2)($_)"});
                    }
                }
            }
        }
    }
    my $model;
    if(exists $run->{"Model"}{1} and $run->{"Model"}{1} eq "Recessive"){
        $model="Recessive";
    }elsif(exists $run->{"Model"}{1} and $run->{"Model"}{1} eq "AutoComHet"){
        $model="Recessive";
    }elsif(exists $run->{"Model"}{1} and $run->{"Model"}{1} eq "Dominance"){
        $model="Dominance";
    }else{
        $model="Dominance";
        print "[WARN] Auto Model: $model\n";
    }
    
    
    open OUT,">$significant_gene" if($DoSigGene eq "TRUE");
    foreach my $Gene(keys %data)
    {
        my ($s0,$s1,$s2)=(0,0,0);
        my ($case_s0,$case_s1,$case_s2)=(0,0,0);
        my ($control_s0,$control_s1,$control_s2)=(0,0,0);
        my ($S1,$S2)=(" "," ");
        foreach my $sample(keys %{$data{$Gene}{"Sample"}})
        {
            my $mut=$data{$Gene}{"Sample"}{$sample};
            my $mark="";
            if($casehash{$sample}){
                $mark="case";
            }
            if($controlhash{$sample}){
                $mark="control";
            }
            if($mut==0){
                $s0++;
                $case_s0++ if($mark eq "case");
                $control_s0++ if($mark eq "control");
            }
            
            if($mut==1){
                $s1++;
                $S1.="$sample,";
                $case_s1++ if($mark eq "case");
                $control_s1++ if($mark eq "control");
            }
            if($mut==2){
                $s2++;
                $S2.="$sample,";
                $case_s2++ if($mark eq "case");
                $control_s2++ if($mark eq "control");
            }
        }
        
        my ($case_freq, $control_freq) = (0, 0);
        my ($p1, $p2, $p3, $p4) = (0, 0, 0, 0);
        my $ControlNum = (exists $run->{"Control_Num_Limit"}{1}) ? $run->{"Control_Num_Limit"}{1} : 50;
        if(@control >= $ControlNum){
            ($p1,$p2,$p3,$p4)=package::p1p2p3p4::run($case_s0,$case_s1,$case_s2,$control_s0,$control_s1,$control_s2);
            $case_freq = (($case_s0+$case_s1+$case_s2) == 0) ? 0 : ($case_s1+$case_s2)/($case_s0+$case_s1+$case_s2);
            $control_freq = (($control_s0+$control_s1+$control_s2) == 0) ? 0 : ($control_s1+$control_s2)/($control_s0+$control_s1+$control_s2);
        }else{
            my $controlDB="";
            $controlDB = $data{$Gene}{"GeneskyExonDB Mutation(0|1|2)"} if(exists $data{$Gene}{"GeneskyExonDB Mutation(0|1|2)"});
            my ($cdb0,$cdb1,$cdb2)=(0,0,0);
            ($cdb0,$cdb1,$cdb2)=split /\|/,$controlDB if($controlDB=~/\d/);
            ($p1,$p2,$p3,$p4)=package::p1p2p3p4::run($case_s0,$case_s1,$case_s2,$cdb0,$cdb1,$cdb2);  
            $case_freq = (($case_s0+$case_s1+$case_s2) == 0) ? 0 : ($case_s1+$case_s2)/($case_s0+$case_s1+$case_s2);
            $control_freq = (($cdb0+$cdb1+$cdb2) == 0) ? 0 : ($cdb1+$cdb2)/($cdb0+$cdb1+$cdb2);
        }
        
        if($DoSigGene eq "TRUE"){
            if($case_freq >= 2*$control_freq){
                print OUT $Gene."\n";
                $GENE -> write($row,0,$Gene,$format->{"red"});
            }elsif($case_freq > $control_freq and $p2 < 0.05){
                print OUT $Gene."\n";
                $GENE -> write($row,0,$Gene,$format->{"red"});
            }else{
                $GENE -> write($row,0,$Gene,$format->{"normal"});
            }
        }else{
            $GENE -> write($row,0,$Gene,$format->{"normal"});
        }

        if(exists($control[0])){
         for my $i(1..@title2-1){
            my $v=" ";$v=$data{$Gene}{$title2[$i]} if(exists $data{$Gene}{$title2[$i]});
            $v=keys %{$data{$Gene}{"Site"}} if($title2[$i] eq "SNV count");
            $v=package::priority::Gene($model,$data{$Gene}{"Priority"}) if($title2[$i] eq "First Priority");
            $v=$s1+$s2 if($title2[$i] eq "Sample count");
            $v="$s0|$s1|$s2" if($title2[$i] eq "Mutation(0|1|2)");
            $v="$case_s0|$case_s1|$case_s2" if($title2[$i] eq "Case(0|1|2)");
            $v="$control_s0|$control_s1|$control_s2" if($title2[$i] eq "Control(0|1|2)");
            $v=$p1 if($title2[$i] eq "p1");
            $v=$p2 if($title2[$i] eq "p2");
            $v=$p3 if($title2[$i] eq "p3");
            $v=$p4 if($title2[$i] eq "p4");
            $v=($s1/2+$s2)/($s0+$s1+$s2) if($title2[$i] eq "Alt Allele Freq");
            $v=$S1 if($title2[$i] eq "Mutation 1");
            $v=$S2 if($title2[$i] eq "Mutation 2");
            $GENE->write($row,$i,$v,$format->{"normal"});
            }
        }else{
         for my $i(1..@title-1){
            my $v=" ";$v=$data{$Gene}{$title[$i]} if(exists $data{$Gene}{$title[$i]});
            $v=keys %{$data{$Gene}{"Site"}} if($title[$i] eq "SNV count");
            $v=package::priority::Gene($model,$data{$Gene}{"Priority"}) if($title[$i] eq "First Priority");
            $v=$s1+$s2 if($title[$i] eq "Sample count");
            $v="$s0|$s1|$s2" if($title[$i] eq "Mutation(0|1|2)");
            $v=($s1/2+$s2)/($s0+$s1+$s2) if($title[$i] eq "Alt Allele Freq");
            $v=$p1 if($title[$i] eq "p1");
            $v=$p2 if($title[$i] eq "p2");
            $v=$p3 if($title[$i] eq "p3");
            $v=$p4 if($title[$i] eq "p4");
            $v=$S1 if($title[$i] eq "Mutation 1");
            $v=$S2 if($title[$i] eq "Mutation 2");
            $GENE->write($row,$i,$v,$format->{"normal"});
            }
        }
        $row++;
    }
    close OUT if($DoSigGene eq "TRUE");
}

1