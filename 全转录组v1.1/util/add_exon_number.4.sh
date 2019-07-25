grep exon $1 | sort -k2,2n -k3,3n > $2/tmp.gtf
perl /home/xudl/perl_scripts/RNA/adapter/gtf_fmt.pl $2/tmp.gtf > $2/tmp.gtf.fmt
for i in A3SS A5SS RI SE;do 
	perl /home/xudl/perl_scripts/RNA/adapter/alter_fmt.pl $2/tmp.gtf.fmt $2/$i.MATS.JC.txt > $3/$i.MATS.JunctionCountOnly.xls;
	perl /home/xudl/perl_scripts/RNA/adapter/alter_fmt.pl $2/tmp.gtf.fmt $2/$i.MATS.JCEC.txt > $3/$i.MATS.ReadsOnTargetAndJunctionCounts.xls;
done;
perl /home/xudl/perl_scripts/RNA/adapter/MXE_alter_fmt.pl $2/tmp.gtf.fmt $2/MXE.MATS.JC.txt > $3/MXE.MATS.JunctionCountOnly.xls;
perl /home/xudl/perl_scripts/RNA/adapter/MXE_alter_fmt.pl $2/tmp.gtf.fmt $2/MXE.MATS.JCEC.txt > $3/MXE.MATS.ReadsOnTargetAndJunctionCounts.xls; 
