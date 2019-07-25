package html_report;

use Cwd 'abs_path';
# use utf8;
use strict;
use Parallel::ForkManager;
use Getopt::Long;
use File::Spec;
use GD;
sub run
{
	my $metadata          = shift;
	my $base              = shift;
	my $config			  = shift;
	$config = abs_path($config);
	my $html_module_dir   = qq{$base->{scriptdir}/html};
	my $Reportpath		  = qq{$metadata->{result}/Report};

	print "======== 请输入用于制作html报告的分组文件名 ========\n";
	print "======== 例如：group.xls则输入group ========\n";
	my $groupfile_id = <STDIN>;
	$groupfile_id =~ s/[\r\n]//g;

	if (-d "$Reportpath/$groupfile_id") {
		Format($Reportpath, $groupfile_id, $html_module_dir, $metadata, $config);

		Html($Reportpath, $groupfile_id, $html_module_dir, $metadata, $config);

	}else{

		die "Can not find file $Reportpath/$groupfile_id, please check it!\n";
	}

	print "html报告已生成！\n";
}

###################################################################子程序

sub Format {
	my $Reportpath   = shift;
	my $groupfile_id = shift;
	my $dir          = shift;
	my $metadata     = shift;
	my $config		 = shift;

	my $intm_result  = qq{$metadata->{intm_result}/html};
	system qq{perl $dir/txt/txt.qc.pl $Reportpath/Statistics $intm_result};
	my @type = ("Relative_Quantitation", "Absolute_Quantitation");

	for (my $i = 0; $i < @type; $i++) {

		if ($type[$i] eq "Relative_Quantitation") {
			system qq{perl $dir/txt/txt.otu.pl $Reportpath/$groupfile_id/$type[$i]/OTU $intm_result $type[$i]};

		}else{

			system qq{perl $dir/txt/txt.otu.IR.pl $Reportpath/$groupfile_id/$type[$i]/OTU $intm_result $type[$i]};
			system qq{perl $dir/txt/txt.Community.IR.pl $Reportpath/$groupfile_id/$type[$i]/Community $intm_result $type[$i]};

		}

		system qq{perl $dir/txt/txt.alpha.pl $Reportpath/$groupfile_id/$type[$i]/AlphaDiversity $intm_result $type[$i]};

		system qq{perl $dir/txt/txt.beta.pl $Reportpath/$groupfile_id/$type[$i]/BetaDiversity $intm_result $type[$i]};

		system qq{perl $dir/txt/txt.picrust.pl $Reportpath/$groupfile_id/$type[$i]/PICRUSt $intm_result $type[$i]};

	}

	system qq{cp $dir/txt/3_5_2.txt $intm_result/txt};
	system qq{cp $dir/txt/info.txt $intm_result/txt};

	
	my $name       = `grep "项目名称" $config |cut -d "=" -f 2`;
	my $project_id = `grep "项目编号" $config |cut -d "=" -f 2`;
	my $work       = `grep "客户单位" $config |cut -d "=" -f 2`;
	my $sample     = `grep "样本来源" $config |cut -d "=" -f 2`;
	my $database   = `grep "参考数据库" $config |cut -d "=" -f 2`;
	my $time 	   = `date "+%Y年%m月%d日"`;
	system qq{echo "$name$project_id$work$sample$database$time" > "$intm_result/txt/info.txt"};

	system qq{cp -r "$dir/template" $intm_result};
	my @tmp = `ls $intm_result/template | grep "txt"`;

	for (my $j = 0; $j < @tmp; $j++) {

		chomp $tmp[$j];
		system qq{sed -i "s!<tab>Relative_Quantitation!<tab>$groupfile_id\/Relative_Quantitation!g" "$intm_result/template/$tmp[$j]"};
		system qq{sed -i "s!<img>Relative_Quantitation!<img>$groupfile_id\/Relative_Quantitation!g" "$intm_result/template/$tmp[$j]"};

		system qq{sed -i "s!<tab>Absolute_Quantitation!<tab>$groupfile_id\/Absolute_Quantitation!g" "$intm_result/template/$tmp[$j]"};
		system qq{sed -i "s!<img>Absolute_Quantitation!<img>$groupfile_id\/Absolute_Quantitation!g" "$intm_result/template/$tmp[$j]"};

		system qq{sed -i "s!<tab>Combine_Analysis!<tab>$groupfile_id\/Combine_Analysis!g" "$intm_result/template/$tmp[$j]"};
		system qq{sed -i "s!<img>Combine_Analysis!<img>$groupfile_id\/Combine_Analysis!g" "$intm_result/template/$tmp[$j]"};

		# system qq{sed -i "s!<tab>txt!<tab>$intm_result\/txt!g" "$intm_result/template/$tmp[$j]"};
	}


	my $level = qq{$metadata->{level}};
	if ($level eq "DNA") {
		system qq{sed -i '16,22d' "$intm_result/template/4_2_1_2.txt"};
		system qq{sed -i '8d' "$intm_result/template/4_2_1_3.txt"};
		system qq{sed -i '7,10d' "$intm_result/template/4_2_5_1.txt"};
		system qq{echo -e "<strong><font color="red" size="3">!!后续多样性分析基于DNA水平拷贝数：Absolute_Quantitation/OTU/otu_copies_unit_DNA.xls。</font></strong>\n<t1>&nbsp;</t1>" >>$intm_result/template/4_2_1_2.txt};	

	}else{
		system qq{sed -i '15,16d' "$intm_result/template/4_2_1_2.txt"};
		system qq{sed -i '9d' "$intm_result/template/4_2_1_3.txt"};
		system qq{sed -i '11,13d' "$intm_result/template/4_2_5_1.txt"};
		system qq{echo -e "<strong><font color="red" size="3">!!后续多样性分析基于样本水平拷贝数：Absolute_Quantitation/OTU/otu_copies_unit_sample.xls。</font></strong>\n<t1>&nbsp;</t1>" >>$intm_result/template/4_2_1_2.txt};
	}


	system qq{sed -i '13,17d' "$intm_result/template/4_4_4_2.txt"} if not -e "$Reportpath/$groupfile_id/Combine_Analysis/BetaDiversity/Lefse/both_but_diff_group_barplot.pdf";

}

sub Html {

	my $Reportpath   = shift;
	my $groupfile_id = shift;
	my $dir          = shift;
	my $metadata     = shift;
	my $config		 = shift;

	my $out_dir      = $Reportpath;
	my $intm_result  = qq{$metadata->{intm_result}/html};
	my $tem_dir		 = "$intm_result/template";
	my $catalog      = "$tem_dir/catalog.txt";
	my $project_id   = `grep "项目编号" $config |cut -d "=" -f 2`;
	my $name       	 = `grep "项目名称" $config |cut -d "=" -f 2`;
	my $time         = `date "+%Y%m%d"`;
	$name =~ s/[\r\n]//g;
	$project_id =~ s/[\r\n]//g;
	$time =~ s/[\r\n]//g;
	my $index_html	 = qq{$out_dir/$project_id\_$name\_结题报告_$time\.html};
	#复制css、js、img等文件
	system "cp -r $dir/html/ $out_dir/";
	system "cp -r $tem_dir/IMG $out_dir/html/" if(-e "$tem_dir/IMG");

	#填充进index.html
	create_index($catalog, $index_html, $out_dir, $tem_dir, $groupfile_id, $intm_result);
}


###################################################################### Html调用的子函数

sub create_index{
	my ($catalog, $index, $out_dir, $tem_dir, $group_id, $intm_result) = @_;

	my $head = '<!DOCTYPE html>
	<html lang="en">
	<head>
		<meta charset="UTF-8">
		<title>项目报告</title>
		<link rel="stylesheet" href="./html/css/bootstrap-theme.css">
		<link rel="stylesheet" href="./html/css/bootstrap.min.css">
		<link rel="stylesheet" href="./html/css/initialize.css">
		<script src="./html/js/jquery-3.2.1.min.js"></script>
		<script src="./html/js/bootstrap.js"></script>
	</head>
	<body>';

	my $tail = '<a href="#top">
		<div class="return_top">
			<span class="glyphicon glyphicon-arrow-up" aria-hidden="true"></span>
		</div>
	</a>
	<script src="./html/js/action.js"></script>
	</body>
	</html>';

	my $top = '
	<div class="top_fixed">
		<div class="top_bk" id="top">
		    <div class="top_nav">
		        <img src="./html/IMG/logo.png" alt="genesky">
		        <div class="top_right">
		            <div>
		                <img src="html/IMG/tel.png">
		                <span>咨询热线：<a>400-065-6886</a></span>
		                <img src="html/IMG/url.png">
		                <span>网址：<a href="http://www.geneskybiotech.com" target="_blank" style="color:#337ab7;">http://www.geneskybiotech.com</a></span>
		            </div>
		            <div>
		                <img src="html/IMG/address.png">
		                <span>地址：<a>上海市浦东新区康桥路787号9号楼</a></span>
		            </div>
		        </div>
		    </div>
		</div>
		<div class="top_catalog">
		    <ul>';
	my $down = '<div class="body_box">
			<div class="main_body">';
	
	###---------3_5_2.txt-----------
	system qq{mkdir -p $intm_result/txt} if not -d qq{$intm_result/txt};
	my $in = "$tem_dir/TXT/3_5_2.txt";
	my $out = ">$intm_result/txt/3_5_2.txt";
	my $hash = read_3_5_2_result($out_dir, $group_id);
	open FILE, $in;
	open OUT, $out;
	while (<FILE>){

		chomp;
		if (/\t/){
			my ($key,$value) = split /\t/;
			if (exists $hash->{$key}){
				print OUT qq{$key\t√\n};
			}else{
				print OUT qq{$key\t \n};
			}
		}else{
			print OUT qq{$_\n};
		}
	}
	close OUT;
	close FILE;
	

	##----------catalog.txt----------------
	system qq{mkdir -p $intm_result/txt} if not -d qq{$intm_result/txt};
	open CATA, $catalog;
	open CATALOG,">$intm_result/txt/catalog_tmp.txt";
	my $hash2 = read_catalog_result($out_dir, $group_id);
	while (<CATA>){
		chomp;
		my $key = (split /\t/, $_)[1];                 
		next if (not -exists $hash2->{$key});	
		print CATALOG qq{$_\n};
	}
	close CATA;
	close CATALOG;
	change_catalog_stepid("$intm_result/txt/catalog_tmp.txt", "$intm_result/txt/catalog.txt");


	open CATA,"$intm_result/txt/catalog.txt";
	my $line2 = <CATA>;

	my $tmp = '';
	while (my $line1 = $line2){

		$line2 = <CATA>;
		$line1 =~s /[\r\n]//g;
		$line2 =~s /[\r\n]//g if (defined $line2);
		
		
		my($id, $first, $content, $file) = split "\t", $line1 ;
		my $chapter = (split "_", $first)[0];
		my $first_ = $first;
		$first =~s /_/./g;
		my $title = "$first.$content";

		my $whole = read_file("$tem_dir/$file", $out_dir, $tem_dir, $intm_result) if (defined $file);

		if($id eq "S1"){

			if(defined $file){
				$top  .= '<li><a href="#">'.$content.'</a></li>';
				$down .= '<div class="part-'.$first_.' content_bk"><div class="content-'.$chapter.'">'.$whole.'</div></div>';
			}else{
				$top  .= '<li><a href="#" class="line_through">'.$content.'</a></li>';
			}

	
		}elsif($id eq "S2"){

			$top.='<li><a href="#">'.$content.'</a></li>';
			$down.='<div class="part-'.$first_.'">
				<div class="body_left_box">
					<div class="catalog-'.$first_.' body_left">
						<h4 class="catalog-title"><span class="glyphicon glyphicon-menu-down" aria-hidden="true"></span>&nbsp;&nbsp;&nbsp;'.$title.'</h4>
							<ul>';


		}elsif($id eq "S3"){

			if(defined $file){

				$tmp .= '<div class = "content-'.$first_.'" id="content-'.$first_.'">'.$whole.'</div>';

				$down.='<li class="catalog-'.$first_.'"><a href="#content-'.$first_.'"><span class="glyphicon glyphicon-menu-right" aria-hidden="true"></span>&nbsp;&nbsp;&nbsp;'.$title.'</a></li>';
			}else{

				$down.='<li class="catalog-'.$first_.'"><a href="#" class="line_through"><span class="glyphicon glyphicon-menu-right" aria-hidden="true"></span>&nbsp;&nbsp;&nbsp;'.$title.'</a></li>';
			}


		}elsif($id eq "S4"){

			if(defined $file){
				$tmp .= '<div class = "content-'.$first_.'" id = "content-'.$first_.'">'.$whole;
				$down.='<li class="catalog-'.$first_.'"><a href="#content-'.$first_.'"><span class="glyphicon glyphicon-menu-right" aria-hidden="true"></span>&nbsp;&nbsp;&nbsp;'.$title.'</a><ul class="catalog_down">';
			}else{
				$down.='<li class="catalog-'.$first_.'"><a href="#" class="line_through"><span class="glyphicon glyphicon-menu-right" aria-hidden="true"></span>&nbsp;&nbsp;&nbsp;'.$title.'</a><ul class="catalog_down">';
			}

		}elsif($id eq "S5"){

			$tmp .= '<div class = "content-'.$first_.'" id = "content-'.$first_.'">'.$whole.'</div>' if (defined $file);
			$down.='<li class="catalog-'.$first_.'"><a href="#content-'.$first_.'">'.$title.'</a></li>';


		}elsif($id eq "S6"){

			$tmp .= '<div class = "content-'.$first_.'" id = "content-'.$first_.'">'.$whole.'</div></div>' if (defined $file);
			$down.='<li class="catalog-'.$first_.'"><a href="#content-'.$first_.'">'.$title.'</a></li></ul></li>';
			

		}elsif($id eq "S7"){

			$tmp .= '<div class = "content-'.$first_.'" id = "content-'.$first_.'">'.$whole.'</div></div>' if (defined $file);
			$down.='<li class="catalog-'.$first_.'"><a href="#content-'.$first_.'">'.$title.'</a></li></ul></li></ul></div></div><div class="body_right_box"><div class="content-'.$chapter.'">'.$tmp.'</div></div></div>';
			$tmp = '';


		}elsif($id eq "S8"){

			if(defined $file){
				$tmp .= '<div class = "content-'.$first_.'" id = "content-'.$first_.'">'.$whole.'</div>';
				$down.='<li class="catalog-'.$first_.'"><a href="#content-'.$first_.'"><span class="glyphicon glyphicon-menu-right" aria-hidden="true"></span>&nbsp;&nbsp;&nbsp;'.$title.'</a></li>
								</ul>
							</div>
						</div>
						<div class="body_right_box"><div class="content-'.$chapter.'">'.$tmp.'
						</div></div></div>
						';
				$tmp = '';
			}else{
				$down.='<li class="catalog-'.$first_.'"><a href="#" class="line_through"><span class="glyphicon glyphicon-menu-right" aria-hidden="true"></span>&nbsp;&nbsp;&nbsp;'.$title.'</a></li>
							</ul>
						</div>
					</div>
					<div class="body_right_box"><div class="content-'.$chapter.'">'.$tmp.'</div></div>
				</div>';
				$tmp = '';
			}

		
		}elsif($id eq "S9"){  

			if(defined $file){
				$tmp .= '<div class = "content-'.$first_.'" id = "content-'.$first_.'">'.$whole.'</div>';
				$down.='<li class="catalog-'.$first_.'"><a href="#content-'.$first_.'"><span class="glyphicon glyphicon-menu-right" aria-hidden="true"></span>&nbsp;&nbsp;&nbsp;'.$title.'</a><ul class="catalog_down">';
			}else{
				$down.='<li class="catalog-'.$first_.'"><a href="#" class="line_through"><span class="glyphicon glyphicon-menu-right" aria-hidden="true"></span>&nbsp;&nbsp;&nbsp;'.$title.'</a><ul class="catalog_down">';
			}


		}elsif($id eq "S10"){   

			$tmp .= '<div class = "content-'.$first_.'" id = "content-'.$first_.'">'.$whole.'</div></div>' if (defined $file);
			$down.='<li class="catalog-'.$first_.'"><a href="#content-'.$first_.'">'.$title.'</a></li></ul></li>';


		}elsif($id eq "S12"){   

			$tmp .= '<div class = "content-'.$first_.'" id = "content-'.$first_.'">'.$whole.'</div></div></div>' if (defined $file);
			$down.='<li class="catalog-'.$first_.'"><a href="#content-'.$first_.'">'.$title.'</a></li></ul></li>';


		}elsif($id eq "S13"){  

			$tmp .= '<div class = "content-'.$first_.'" id = "content-'.$first_.'">'.$whole.'</div></div>' if (defined $file);
			$down.='<li class="catalog-'.$first_.'"><a href="#content-'.$first_.'">'.$title.'</a></li></ul></li></ul></li>';


		}elsif($id eq "S14"){  
			if(defined $file){
				$tmp .= '<div class = "content-'.$first_.'" id = "content-'.$first_.'">'.$whole.'</div>';
				$down.='<li class="catalog-'.$first_.'"><a href="#content-'.$first_.'"><span class="glyphicon glyphicon-menu-right" aria-hidden="true"></span>&nbsp;&nbsp;&nbsp;'.$title.'</a></li></ul></li></ul></li>
								</ul>
							</div>
						</div>
						<div class="body_right_box"><div class="content-'.$chapter.'">'.$tmp.'
						</div></div></div></div>
						';
				$tmp = '';
			}else{
				$down.='<li class="catalog-'.$first_.'"><a href="#" class="line_through"><span class="glyphicon glyphicon-menu-right" aria-hidden="true"></span>&nbsp;&nbsp;&nbsp;'.$title.'</a></li>
							</ul>
						</div>
					</div>
					<div class="body_right_box"><div class="content-'.$chapter.'">'.$tmp.'</div></div>
				</div>';
				$tmp = '';
			}

		}

	}
	

	$top .= '</ul>
			</div>
	</div>
	';
	$down .= '</div>
	</div>
	';	
	close CATA;

	open INDEX, ">$index";
	print INDEX $head, $top, $down, $tail;
	close INDEX;

}


sub read_file{
	my ($file, $out_dir, $tem_dir, $intm_result) = @_;
	my $last = '';

	open TEM, $file or die " -- [ERR]不能读取 $file \n";
	while (my $line = <TEM>){ 
		$line =~s/[\r\n]//g;
		if($line =~ /^<t1>(.*)<\/t1>/){				
			$line =~ s/<t1>/<p class="paragraph">/;
			$line =~ s/<\/t1>/<\/p>/;
		}
		if($line =~ /^<t2>(.*)<\/t2>/){				
			$line =~ s/<t2>/<p class="text_center"; style= "font-size:15px">/;
			$line =~ s/<\/t2>/<\/p>/;
		}
		if($line =~ /^<t3>(.*)<\/t3>/){				
			$line =~ s/<t3>/<p class="text_left"; style="font-family:楷体; font-size:15px">/;
			$line =~ s/<\/t3>/<\/p>/;
		}
		if($line =~ /^<t4>(.*)<\/t4>/){				
			$line =~ s/<t4>/<p class="text_left"; style="background:Gainsboro">/;
			$line =~ s/<\/t4>/<\/p>/;
		}
		if($line =~ /<tab(.*)>(.+)<\/tab>$/){
			my $align = $1;
			my $tmp = $2;

			my $dir;
			if ($tmp =~ /^#/) {
				$dir = $tem_dir;
			}elsif($tmp =~ /^txt/){
				$dir = $intm_result;
			}else{
				$dir = $out_dir;
			}

			$tmp =~s /^#//;
			# print "$tmp\n";

			if (`ls $dir/$tmp | wc -l` > 0){
				if ($tmp =~ /index.txt/ and $dir eq $tem_dir){

					# print "-----------1------------\n";
					$line = create_table_fengmian("$dir/$tmp", $intm_result);

				}elsif($tmp =~ /^txt/){
					# print "-----------2------------\n";
					$line = create_table_normal($dir, $tmp, $align);
				}else{
					# print "$tmp\n";
					# print "-----------3------------\n";
					$line = create_table_normal($dir, $tmp, $align);
				}
			}
		}
		if($line =~ /<img>(.*)<\/img>$/){
			my $file_dir = $1;
			my $dir = ($file_dir =~ /^#/)? $tem_dir: $out_dir;
			my $prefix = ($file_dir =~ /^#/)? "html/" : "";
			$file_dir =~ s/^#//;
			
			my $filedir = (File::Spec->splitpath(File::Spec->rel2abs("$dir/$file_dir")))[1];
			my $outdir  = (File::Spec->splitpath(File::Spec->rel2abs("$dir/*")))[1];
			my $filename = (File::Spec->splitpath(File::Spec->rel2abs("$dir/$file_dir")))[2];
			$filedir =~ s/\*\///g;
			my $num = `find $filedir -name $filename | wc -l`;

			if ( $num > 0 ){
				if ($num > 1){

					# 获取图片
					my @files_dir = `find $filedir -name $filename` ;
					die " -- [ERR] $file 文件设置错误:无满足的 $file_dir\n" unless (@files_dir>0);
					# 样本选择
					$line = '<form class="form-inline"><label class="select_lab">样本选择：</label><select class="select_category input-sm form-control">';

					my %files;
					foreach (@files_dir){
						chomp $_;
						$_ =~s /$outdir//;
						$files{$_} =  $_;
						
					}

					my $first;
					foreach (sort keys %files){
						$first++;
						$line .='<option value="genus_'.$first.'">'.$_.'</option>';
					}
					$line .= '</select>';

					$first = 0;
					if ($file_dir =~ /\.png|\.jpg|\.gif/){
						foreach (sort keys %files){
							my $image = GD::Image->new("$dir/$files{$_}");
							my $width = $image->width;
							$first++;
							if ($width > 800){
								$line .= ($first == 1)? '<img class="img_emb genus_1" style="display:block;width:800px;" src="'."$prefix$files{$_}".'" />' : '<img class="img_emb genus_'.$first.'"  style="width:800px;" src="'."$prefix$files{$_}".'"/>';
							}else{
								$line .= ($first == 1)? '<img class="img_emb genus_1" style="display:block;" src="'."$prefix$files{$_}".'" />': '<img class="img_emb genus_'.$first.'" src="'."$prefix$files{$_}".'" />';
							}
						}
						$line .= '</form>';	
					}elsif($file_dir =~ /\.pdf|\.html|\.htm/){
						my $tmp = ($file_dir =~ /\.pdf/)? "embed":"iframe";
						foreach (sort keys %files){
							$first++;
							if ($tmp =~ /embed/){
								$line .= ($first == 1)?'<'.$tmp.' class="img_emb genus_1" style="display:block;" src="'."$prefix$files{$_}".'" type="application/pdf" width="800" height="800" >':'<'.$tmp.' class="img_emb genus_'.$first.'" src="'."$prefix$files{$_}".'" type="application/pdf" width="800" height="800">';
							}else{
								$line .= ($first == 1)?'<'.$tmp.' class="img_emb genus_1" style="display:block;" src="'."$prefix$files{$_}".'"  width="800" height="800" scroll="yes" ></iframe>':'<'.$tmp.' class="img_emb genus_'.$first.'" src="'."$prefix$files{$_}".'"  width="800" height="800" scroll="yes"></iframe>';
							}
							
						}
						$line .= '</form>';
					}	
				}else{
					if ($file_dir =~ /\.png|\.jpg|\.gif/){
						die " -- [ERR] 无此图片 $dir/$file_dir\n" if (not -e "$dir/$file_dir");
						my $image = GD::Image->new("$dir/$file_dir");
						my $width = $image->width;
						$line = ($width > 800)? '<img src="'."$prefix$file_dir".'" style="width:800px;"/>' : '<img src="'."$prefix$file_dir".'"/>';
					}elsif ($file_dir =~ /\.pdf/){
						# 获取图片
						$file_dir = `find $filedir -name $filename` ;
						$file_dir =~s /$outdir//;
						$line = '<embed src="'."$prefix$file_dir".'" width="800" height="800" />';
					}elsif ($file_dir =~ /\.html|\.htm/){
						$line = '<iframe src="'."$prefix$file_dir".'" width="800" height="800" scroll="yes"></iframe>';
					}
				}
			}
		}
		$last .= $line;
	}
	close TEM;
	
	return $last;
}


sub create_table_fengmian{
	my ($file, $intm_result) = @_;
	#读取项目信息
	my $info = "$intm_result/txt/info.txt";
	open INFO, $info or die " -- [ERR] 不能读取 $info \n";
	my @infos;
	while (<INFO>){
		next if ($_ =~ /#/);
		$_ =~s/[\r\n]//g;
		next if ($_ !~ /\S+/);
		$_ =~s/(.*):(.*)/$2/;
		push @infos, $_;
	}
	close INFO;
	#封面表格输出
	my $table;
	open FILE, $file  or die " -- [ERR] 不存在$file\n";
	my $row;
	while (<FILE>){
		next if ($_ =~ /^#/);
		$_ =~s /[\r\n]//g;
		$row++;
		$table .= ($row == 1)? '<tbody><tr>' : '<tr>' ;
		my $count = 0; #单元格
		if ($_ =~/^\[(.*)\]$/ ){#列表
			$_ =~s /\]|\[//g;
			my @array = split "\t", $_;
			foreach my $td(@array){
				$count ++;
				if ($td =~ /<img>(.*)<\/img>$/){
					my $file_dir = $1;
					$file_dir =~s /^#//;
					my ($front) = $td =~ /(.*)<img>/;
					$td = $front.'<img src="./html/'.$file_dir.' "alt></img>'
				};
				$table .= '<td colspan ="4" ><ul id = "project-msg"><li>'.$td.'</li>' if ($count == 1);
				$table .= '<li>'.$td.'</li>' if ($count !=1 and $count !=3);
				if ($count == 3){
					$td = (@infos > 0)? shift @infos: $td;
					$table .= '<li>'.$td.'</li>' ;
				}
			}
			$table .= '</ul></td>';
		}else{	
			my @array = split "\t", $_;
			my $num = scalar @array;
			foreach my $td(@array){
				my $rowspan = $td =~s /\*//g;
				$rowspan = ($rowspan eq '')? 1: $rowspan;
				if ($num == 1){
					$table .= '<td class="text-center" colspan="4" rowspan="'.$rowspan.'">'.$td.'</td>' ;
				}elsif ($num == 2){
					$count++;
					if ($count == 1){
						$table .= '<td class="td-back" rowspan="'.$rowspan.'">'.$td.'</td>' ;
					}else{
						if (@infos > 1){#因为最后一项是 时间
							my $tmp = shift @infos;
							$td = ($tmp =~ /\S+/)? $tmp : $td;
						}
						$table .= '<td class="text-center" colspan="3" rowspan="'.$rowspan.'">'.$td.'</td>'
					 }
				}else{
					$count++;
					$table .= '<td rowspan="'.$rowspan.'">'.$td.'</td>' if ($count == 1 or ($count%4) != 1);
					$table .= '<tr><td rowspan="'.$rowspan.'">'.$td.'</td>' if ($count != 1 and ($count%4)== 1);
					$table .= '</tr>' if (($count%4) == 0);
					$table .= '</tr>' if (($count%4) != 0 and  $count == $num );
				}
			}
		}
		$table .= '</tr>';
		
	}
	my $head = '<div class="table-div">
    <table class="special-table table-all">
	';
	my $tail = '</tbody></table></div>';
	$table = "$head$table$tail";
	return $table;
}


sub create_table_normal{
	my ($dir, $tmp, $align) = @_;
	my $file = "$dir/$tmp";
	my $table = '';
	$align = ($align =~ /left/)? "left": "center";
	if ($file =~ /\*/){
		my $same = '';
		foreach(split /\/+/, $tmp){
			$same .= "$_/" if ($_ !~ /\*/);
			last if ($_ =~ /\*/);
		}
		# 获取表格
		my @files_dir = `ls $file` ;
		die " -- [ERR] 文件设置错误:无满足的$file\n" unless (@files_dir>0);
		# 获取表格列数
		my %cols = count_col($files_dir[0]);
		# 获取表格内容
		my %files;
		my $content = '';
		my $first;
		foreach my $file(@files_dir){
			$first++;
			chomp $file;
			$content  .= table($file, $align, \%cols, $first);
			$file =~s /$dir\/*$same\/*(.*)/$1/;
			$files{$file} = ($same eq '')? $file : "$same/$file";
		}
		
		$table = '<form class="form-inline"><label class="select_lab">样本选择：</label><select class="select_category input-sm form-control">';
		$first = 0;
		foreach (sort keys %files){
			$first++;
			$table .='<option value="genus_'.$first.'">'.$_.'</option>';
		}
		$table .= '</select>';
		$table .= $content;
		$table .= '</form>';
	}else{
		#获取表格列数
		my %cols = count_col($file);
		my $max = (sort {$b <=> $a} keys %cols)[0];
		# 获取表格内容
		$table = table($file, $align, \%cols);	
	}
	return $table;
}


sub table{
	my ($file, $align, $cols, $first) = @_;
	# print "$file\n";
	my $tmp;
	if (not defined $first){

		my @tmpf = split /\//, $file; 
		if ($tmpf[@tmpf-1] =~ /3_5_2\.txt$/ or ($tmpf[@tmpf-1] =~ /7_2\.txt$/)) {
			$tmp = '<div class="scroll-div" style="overflow:auto;height:600x;"><table class="scroll-table" ">';

		}elsif($tmpf[@tmpf-1] =~ /txt$/){
			
			$tmp = '<div class="scroll-div" style="overflow:auto;height:260px;"><table class="scroll-table" ">';

		}else{
			$tmp = '<div class="scroll-div" style="overflow:auto;height:300px;"><table class="scroll-table" ">';
		}
		

	}elsif($first == 1){
		$tmp ='<div class="img_emb genus_'.$first.'" style="display:block;overflow:auto;height:300px;"><table class="scroll-table" ">';

	}else{
		$tmp = '<div class="img_emb genus_'.$first.'" style="overflow:auto;height:300px;"><table class="scroll-table" ">';#img_emb scroll-div_emb

	}	
	# 获取表格内容
	my $head = ($file =~ /\/+2.txt$/)?'<div class="table-div"><table class="style-table" ">': $tmp; #style ="width:100%;
	my $tail = '</tbody></table></div>';
	my ($table, $row);
	my $max = (sort {$b <=> $a} keys %{$cols})[0];
	open FILE, $file;
	while (<FILE>){
		next if ($_ =~ /^#/);
		chomp $_;
		$_ =~s /[\r\n]//;
		$row++;
		my @tds = split /\t/, $_;
		if (keys %{$cols} == 1){#规整的 标题和内容
			my $count = 0;
			$table .= '<thead><tr>' if ($row == 1);
			$table .= '<tbody><tr>' if ($row == 2);
			$table .= '<tr>' if ($row != 2 and $row != 1 );
			foreach (@tds){
				$count++;
				if ($row == 1){
					$table .= ($count == 1 )? '<th style= "text-align:'.$align.'">'.$_.'</th>' : '<th style="text-align:'.$align.'">'.$_.'</th>';
				}else{
					$table .= ($count == 1 )? '<td style= "text-align:'.$align.'">'.$_.'</td>' : '<td style="text-align:'.$align.'">'.$_.'</td>';
				}
			}
			$table .= ($row == 1)?  '</tr></thead>' : '</tr>';

		}else{
			$table .= ($row == 1)? '<tbody><tr>' : '<tr>' ;
			my $count = 0;# td index
			if ($_ =~ /^\[/ and $_=~/\]$/ ){#列表
				$_ =~s /\[|\]//g;
				foreach my $td(@tds){
					$count ++;
					if ($td =~ /<img>(.*)<\/img>$/){
						my $file_dir = $1;
						$file_dir =~s /^#//;
						$td = '<img src="./html/'.$file_dir.' "alt></img>';
					}
					$table .= '<td colspan ="'.$max.'" style = "text-align: '.$align.'"><ul id = "project-msg"><li>'.$td.'</li>' if ($count == 1);
					$table .= '<li>'.$td.'</li>' if ($count != 1);
				}
				$table .= '</ul</td>';
			}else{
				my $num = scalar @tds;
				foreach my $td(@tds){
					my $rowspan = ($td) =~s /\*//g;
					my $colspan = ($td) =~s /\@//g;
					$rowspan = ($rowspan =~/\d/)? $rowspan : 1;
					$colspan = ($colspan =~/\d/)? $colspan : 1;
					$count ++;
					$table .= '<td rowspan="'.$rowspan.'" colspan="'.($colspan).'" style="text-align:'.$align.'">'.$td.'</td>' if ( $count == 1 or ($count%$max) != 1) ;
					$table .= '<tr><td rowspan="'.$rowspan.'" colspan="'.($colspan).'" style="text-align:'.$align.'">'.$td.'</td>' if ( $count != 1 and ($count%$max) == 1) ;
					$table .= '</tr>' if (($count%$max) == 0);
					$table .= '</tr>' if (($count%$max) != 0 and $count == $num);	
				}
			}
			$table .= '</tr>';
		}
	}
	$table = "$head$table$tail";
	return $table;

}

sub count_col{
	my $file = shift @_;
	my %cols;
	my $max_col = 0;
	open FILE, $file  or die " -- [ERR] 不存在$file\n";
	while (<FILE>){
		next if ($_ =~ /^#/);
		$_ =~s /[\r\n]//g;
		my @array = split /\t/, $_;
		my $num = @array;
		$max_col = ($num > $max_col)? $num : $max_col;
		$cols{$num}++;
	}
	close FILE;
	return %cols;
}


sub read_catalog_result{

	my $out_dir1= shift @_;
	my $group_id = shift @_;
	my %hash = ();

	my @type = ("Statistics", "Absolute_Quantitation", "Relative_Quantitation", "Combine_Analysis");

		$hash{'1'}      = qq{1};
		$hash{'2'}      = qq{2};
		$hash{'3'}      = qq{3};
		$hash{'3_1'}    = qq{3_1};
		$hash{'3_2'}    = qq{3_2};
		$hash{'3_3'}    = qq{3_3};
		$hash{'3_4'}    = qq{3_4};
		$hash{'3_5'}    = qq{3_5};
		$hash{'3_5_1'}  = qq{3_5_1};
		$hash{'3_5_2'}  = qq{3_5_2};
		$hash{'4'}      = qq{4};

	for (my $i = 0; $i < @type; $i++) {

		if ($type[$i] eq "Statistics") {
			$hash{"4_1"}                     = qq{4_1} if -d qq{$out_dir1/Statistics};

		}elsif($type[$i] eq "Absolute_Quantitation"){

			my $out_dir = "$out_dir1/$group_id/Absolute_Quantitation";

			$hash{'4_2'}                     = qq{4_2}  if -d qq{$out_dir};
			$hash{'4_2_1'}                   = qq{4_2_1}  if -d qq{$out_dir/OTU};
			$hash{'4_2_1_1'}                 = qq{4_2_1_1}  if -d qq{$out_dir/OTU};
			$hash{'4_2_1_2'}                 = qq{4_2_1_2}  if -d qq{$out_dir/OTU};
			$hash{'4_2_1_3'}                 = qq{4_2_1_3}  if -d qq{$out_dir/OTU};
			$hash{'4_2_1_4'}                 = qq{4_2_1_4}  if -d qq{$out_dir/OTU};
			$hash{'4_2_1_5'}                 = qq{4_2_1_5}  if -d qq{$out_dir/OTU};

			$hash{'4_2_2'}                   = qq{4_2_2}  if -d qq{$out_dir/Community};
			$hash{'4_2_2_1'}                 = qq{4_2_2_1}  if -d qq{$out_dir/Community/Pieplot};
			$hash{'4_2_2_2'}                 = qq{4_2_2_2}  if -d qq{$out_dir/Community/KronaPlot};
			$hash{'4_2_2_3'}                 = qq{4_2_2_3}  if -d qq{$out_dir/Community/TaxonTree};
			$hash{'4_2_2_4'}                 = qq{4_2_2_4}  if -d qq{$out_dir/Community/Community_Structure};
			$hash{'4_2_2_5'}                 = qq{4_2_2_5}  if -d qq{$out_dir/Community/Bubble};
			$hash{'4_2_2_6'}                 = qq{4_2_2_6}  if -d qq{$out_dir/Community/Barplot};
			$hash{'4_2_2_7'}                 = qq{4_2_2_7}  if -d qq{$out_dir/Community/Treebar};
			$hash{'4_2_2_8'}                 = qq{4_2_2_8}  if -d qq{$out_dir/Community/Heatmap};
			$hash{'4_2_2_9'}                 = qq{4_2_2_9}  if -d qq{$out_dir/Community/Wilcoxon};
			$hash{'4_2_2_10'}                = qq{4_2_2_10}  if -d qq{$out_dir/Community/Kruskal_Wallis};
			$hash{'4_2_2_11'}                = qq{4_2_2_11} if -d qq{$out_dir/Community/ANOVA};
			$hash{'4_2_2_12'}                = qq{4_2_2_12} if -d qq{$out_dir/Community/Gene_Copies_Correction};

			$hash{'4_2_3'}                   = qq{4_2_3}  if -d qq{$out_dir/AlphaDiversity};
			$hash{'4_2_3_1'}                 = qq{4_2_3_1}  if -d qq{$out_dir/AlphaDiversity/DiversityIndex};
			$hash{'4_2_3_2'}                 = qq{4_2_3_2}  if -d qq{$out_dir/AlphaDiversity/DiversityIndex/Pairwise_groups};
			$hash{'4_2_3_3'}  			     = qq{4_2_3_3}  if -d qq{$out_dir/AlphaDiversity/Rarefaction};
			$hash{'4_2_3_4'}                 = qq{4_2_3_4}  if -d qq{$out_dir/AlphaDiversity/RankAbundance};
			$hash{'4_2_3_5'}                 = qq{4_2_3_5}  if -d qq{$out_dir/AlphaDiversity/Specaccum};

			$hash{'4_2_4'}                   = qq{4_2_4}  if -d qq{$out_dir/BetaDiversity};
			$hash{'4_2_4_1'}                 = qq{4_2_4_1}  if -d qq{$out_dir/BetaDiversity/Venn};
			$hash{'4_2_4_2'}                 = qq{4_2_4_2}  if -d qq{$out_dir/BetaDiversity/SamplesTree};
			$hash{'4_2_4_3'}      			 = qq{4_2_4_3}  if -d qq{$out_dir/BetaDiversity/PCA};
			$hash{'4_2_4_4'}     			 = qq{4_2_4_4}  if -d qq{$out_dir/BetaDiversity/PCoA};
			$hash{'4_2_4_5'}     			 = qq{4_2_4_5}  if -d qq{$out_dir/BetaDiversity/NMDS};
			$hash{'4_2_4_6'}     			 = qq{4_2_4_6}  if -d qq{$out_dir/BetaDiversity/PLS_DA};
			$hash{'4_2_4_7'}                 = qq{4_2_4_7}  if -d qq{$out_dir/BetaDiversity/ADONIS};
			$hash{'4_2_4_8'}                 = qq{4_2_4_8}  if -d qq{$out_dir/BetaDiversity/Lefse};
			$hash{'4_2_4_9'}                 = qq{4_2_4_9}  if -d qq{$out_dir/BetaDiversity/BetaNTI};

			$hash{'4_2_5'}                   = qq{4_2_5}  if -d qq{$out_dir/PICRUSt};
			$hash{'4_2_5_1'}                 = qq{4_2_5_1}  if -d qq{$out_dir/PICRUSt};
			$hash{'4_2_5_2'}                 = qq{4_2_5_2}  if -d qq{$out_dir/PICRUSt/cog_barplot};
			$hash{'4_2_5_3'}       		     = qq{4_2_5_3}  if -d qq{$out_dir/PICRUSt/pair_group_different_analysis};
			$hash{'4_2_5_4'}     			 = qq{4_2_5_4}  if -d qq{$out_dir/PICRUSt/all_group_different_analysis};

			$hash{'4_2_6'}                   = qq{4_2_6}  if -d qq{$out_dir/EnvironmentFactors};
			$hash{'4_2_6_1'}                 = qq{4_2_6_1}  if -d qq{$out_dir/EnvironmentFactors/MantelTest};
			$hash{'4_2_6_2'}                 = qq{4_2_6_2}  if -d qq{$out_dir/EnvironmentFactors/BioEnv};
			$hash{'4_2_6_3'}  				 = qq{4_2_6_3}  if -d qq{$out_dir/EnvironmentFactors/Correlation};
			$hash{'4_2_6_4'}                 = qq{4_2_6_4}  if -d qq{$out_dir/EnvironmentFactors/RDA_CCA};

		}elsif($type[$i] eq "Relative_Quantitation"){

			my $out_dir = "$out_dir1/$group_id/Relative_Quantitation";

			$hash{'4_3'}                     = qq{4_3}  if -d qq{$out_dir};
			$hash{'4_3_1'}                   = qq{4_3_1}  if -d qq{$out_dir/OTU};
			$hash{'4_3_1_1'}                 = qq{4_3_1_1}  if -d qq{$out_dir/OTU};
			$hash{'4_3_1_2'}                 = qq{4_3_1_2}  if -d qq{$out_dir/OTU};
			$hash{'4_3_1_3'}                 = qq{4_3_1_3}  if -d qq{$out_dir/OTU};
			$hash{'4_3_1_4'}                 = qq{4_3_1_4}  if -d qq{$out_dir/OTU};

			$hash{'4_3_2'}                   = qq{4_3_2}  if -d qq{$out_dir/Community};
			$hash{'4_3_2_1'}                 = qq{4_3_2_1}  if -d qq{$out_dir/Community/Pieplot};
			$hash{'4_3_2_2'}                 = qq{4_3_2_2}  if -d qq{$out_dir/Community/KronaPlot};
			$hash{'4_3_2_3'}                 = qq{4_3_2_3}  if -d qq{$out_dir/Community/TaxonTree};
			$hash{'4_3_2_4'}                 = qq{4_3_2_4}  if -d qq{$out_dir/Community/Community_Structure};
			$hash{'4_3_2_5'}                 = qq{4_3_2_5}  if -d qq{$out_dir/Community/Bubble};
			$hash{'4_3_2_6'}                 = qq{4_3_2_6}  if -d qq{$out_dir/Community/Barplot};
			$hash{'4_3_2_7'}                 = qq{4_3_2_7}  if -d qq{$out_dir/Community/Treebar};
			$hash{'4_3_2_8'}                 = qq{4_3_2_8}  if -d qq{$out_dir/Community/Heatmap};
			$hash{'4_3_2_9'}                 = qq{4_3_2_9}  if -d qq{$out_dir/Community/Wilcoxon};
			$hash{'4_3_2_10'}                = qq{4_3_2_10}  if -d qq{$out_dir/Community/Kruskal_Wallis};
			$hash{'4_3_2_11'}                = qq{4_3_2_11} if -d qq{$out_dir/Community/Metastats};
			$hash{'4_3_2_12'}                = qq{4_3_2_12} if -d qq{$out_dir/Community/ANOVA};


			$hash{'4_3_3'}                   = qq{4_3_3}  if -d qq{$out_dir/AlphaDiversity};
			$hash{'4_3_3_1'}                 = qq{4_3_3_1}  if -d qq{$out_dir/AlphaDiversity/DiversityIndex};
			$hash{'4_3_3_2'}                 = qq{4_3_3_2}  if -d qq{$out_dir/AlphaDiversity/DiversityIndex/Pairwise_groups};
			$hash{'4_3_3_3'}  			     = qq{4_3_3_3}  if -d qq{$out_dir/AlphaDiversity/Rarefaction};
			$hash{'4_3_3_4'}  			     = qq{4_3_3_4}  if -d qq{$out_dir/AlphaDiversity/Shannon};
			$hash{'4_3_3_5'}                 = qq{4_3_3_5}  if -d qq{$out_dir/AlphaDiversity/RankAbundance};
			$hash{'4_3_3_6'}                 = qq{4_3_3_6}  if -d qq{$out_dir/AlphaDiversity/Specaccum};


			$hash{'4_3_4'}                   = qq{4_3_4}  if -d qq{$out_dir/BetaDiversity};
			$hash{'4_3_4_1'}                 = qq{4_3_4_1}  if -d qq{$out_dir/BetaDiversity/Venn};
			$hash{'4_3_4_2'}                 = qq{4_3_4_2}  if -d qq{$out_dir/BetaDiversity/SamplesTree};
			$hash{'4_3_4_3'}      			 = qq{4_3_4_3}  if -d qq{$out_dir/BetaDiversity/PCA};
			$hash{'4_3_4_4'}     			 = qq{4_3_4_4}  if -d qq{$out_dir/BetaDiversity/PCoA};
			$hash{'4_3_4_5'}     			 = qq{4_3_4_5}  if -d qq{$out_dir/BetaDiversity/NMDS};
			$hash{'4_3_4_6'}     			 = qq{4_3_4_6}  if -d qq{$out_dir/BetaDiversity/PLS_DA};
			$hash{'4_3_4_7'}                 = qq{4_3_4_7}  if -d qq{$out_dir/BetaDiversity/ADONIS};
			$hash{'4_3_4_8'}                 = qq{4_3_4_8}  if -d qq{$out_dir/BetaDiversity/Lefse};
			$hash{'4_3_4_9'}                 = qq{4_3_4_9}  if -d qq{$out_dir/BetaDiversity/BetaNTI};

			$hash{'4_3_5'}                   = qq{4_3_5}  if -d qq{$out_dir/PICRUSt};
			$hash{'4_3_5_1'}                 = qq{4_3_5_1}  if -d qq{$out_dir/PICRUSt};
			$hash{'4_3_5_2'}                 = qq{4_3_5_2}  if -d qq{$out_dir/PICRUSt/cog_barplot};
			$hash{'4_3_5_3'}       		     = qq{4_3_5_3}  if -d qq{$out_dir/PICRUSt/pair_group_different_analysis};
			$hash{'4_3_5_4'}     			 = qq{4_3_5_4}  if -d qq{$out_dir/PICRUSt/all_group_different_analysis};

			$hash{'4_3_6'}                   = qq{4_3_6}  if -d qq{$out_dir/EnvironmentFactors};
			$hash{'4_3_6_1'}                 = qq{4_3_6_1}  if -d qq{$out_dir/EnvironmentFactors/MantelTest};
			$hash{'4_3_6_2'}                 = qq{4_3_6_2}  if -d qq{$out_dir/EnvironmentFactors/BioEnv};
			$hash{'4_3_6_3'}                 = qq{4_3_6_3}  if -d qq{$out_dir/EnvironmentFactors/Correlation};
			$hash{'4_3_6_4'}                 = qq{4_3_6_4}  if -d qq{$out_dir/EnvironmentFactors/RDA_CCA};

		}elsif($type[$i] eq "Combine_Analysis"){

			my $out_dir = "$out_dir1/$group_id/Combine_Analysis";

			$hash{'4_4'}                     = qq{4_4}  if -d qq{$out_dir};	
			$hash{'4_4_1'}                   = qq{4_4_1}  if -d qq{$out_dir};	

			$hash{'4_4_2'}                   = qq{4_4_2}  if -d qq{$out_dir/Community};
			$hash{'4_4_2_1'}                 = qq{4_4_2_1}  if -d qq{$out_dir/Community/Barplot};
			$hash{'4_4_2_2'}                 = qq{4_4_2_2}  if -d qq{$out_dir/Community/Wilcoxon};
			$hash{'4_4_2_3'}                 = qq{4_4_2_3}  if -d qq{$out_dir/Community/Kruskal_Wallis};
			$hash{'4_4_2_4'}                 = qq{4_4_2_4}  if -d qq{$out_dir/Community/ANOVA};

			$hash{'4_4_3'}                   = qq{4_4_3}  if -d qq{$out_dir/AlphaDiversity};
			$hash{'4_4_3_1'}                 = qq{4_4_3_1}  if -d qq{$out_dir/AlphaDiversity/DiversityIndex};

			$hash{'4_4_4'}                   = qq{4_4_4}  if -d qq{$out_dir/BetaDiversity};
			$hash{'4_4_4_1'}                 = qq{4_4_4_1}  if -d qq{$out_dir/BetaDiversity/ADONIS};
			$hash{'4_4_4_2'}                 = qq{4_4_4_2}  if -d qq{$out_dir/BetaDiversity/Lefse};

			$hash{'4_4_5'}                   = qq{4_4_5}  if -d qq{$out_dir/PICRUSt};
			$hash{'4_4_5_1'}                 = qq{4_4_5_1}  if -d qq{$out_dir/PICRUSt/pair_group_different_analysis};
			$hash{'4_4_5_2'}                 = qq{4_4_5_2}  if -d qq{$out_dir/PICRUSt/all_group_different_analysis};			

		}
		
	}

	$hash{'5'}      = qq{5};
	$hash{'5_1'}    = qq{5_1};
	$hash{'5_2'}    = qq{5_2};
	$hash{'5_3'}    = qq{5_3};
	$hash{'6'}      = qq{6};
	$hash{'7'}      = qq{7};
	$hash{'7_1'}    = qq{7_1};
	$hash{'7_2'}    = qq{7_2};
	$hash{'8'}      = qq{8};
	$hash{'8_1'}    = qq{8_1};
	$hash{'8_2'}    = qq{8_2};
	$hash{'8_3'}    = qq{8_3};

	return \%hash;
}


sub read_3_5_2_result{

	my $out_dir1= shift @_;
	my $group_id = shift @_;

	my @type = ("Statistics", "Absolute_Quantitation", "Relative_Quantitation", "Combine_Analysis");
	my %hash = ();
	for (my $i = 0; $i < @type; $i++) {

		if ($type[$i] eq "Statistics") {
			$hash{"数据质控与统计（Statistics）"}                 = qq{4_1}    if -d qq{$out_dir1/Statistics};
			$hash{"有效及优化序列统计"}                 		  = qq{4_1}    if -d qq{$out_dir1/Statistics};
		}elsif($type[$i] eq "Absolute_Quantitation"){

			my $out_dir = "$out_dir1/$group_id/Absolute_Quantitation";

			$hash{'绝对定量数据分析（Absolute_Quantitation）'}    = qq{4_2}  if -d qq{$out_dir};
			$hash{'OTU分类分析'}                                  = qq{4_2_1}  if -d qq{$out_dir/OTU};
			$hash{'OTU聚类分析'}                                  = qq{4_2_1_1}  if -d qq{$out_dir/OTU};
			$hash{'绝对拷贝数计算'}                               = qq{4_2_1_2}  if -d qq{$out_dir/OTU};
			$hash{'分类学分析'}                                   = qq{4_2_1_3}  if -d qq{$out_dir/OTU};
			$hash{'物种注释统计'}                                 = qq{4_2_1_4}  if -d qq{$out_dir/OTU};
			$hash{'进化树分析'}                                   = qq{4_2_1_5}  if -d qq{$out_dir/OTU};

			$hash{'群落组成分析（Community）'}                    = qq{4_2_2}  if -d qq{$out_dir/Community};
			$hash{'单样本物种组成饼图（Pieplot）'}                = qq{4_2_2_1}  if -d qq{$out_dir/Community/Pieplot};
			$hash{'单样本多级物种组成图（KronaPlot）'}            = qq{4_2_2_2}  if -d qq{$out_dir/Community/KronaPlot};
			$hash{'分类学系统组成树（TaxonTree）'}                = qq{4_2_2_3}  if -d qq{$out_dir/Community/TaxonTree};
			$hash{'物种总丰度柱状图（Community_Structure）'}      = qq{4_2_2_4}  if -d qq{$out_dir/Community/Community_Structure};
			$hash{'丰度分布Bubble图（Bubble）'}                   = qq{4_2_2_5}  if -d qq{$out_dir/Community/Bubble};
			$hash{'多样本物种组成柱状图（Barplot）'}              = qq{4_2_2_6}  if -d qq{$out_dir/Community/Barplot};
			$hash{'多样本物种组成聚类图（Treebar）'}              = qq{4_2_2_7}  if -d qq{$out_dir/Community/Treebar};
			$hash{'Heatmap热图（Heatmap）'}                       = qq{4_2_2_8}  if -d qq{$out_dir/Community/Heatmap};
			$hash{'Wilcoxon差异分析（Wilcoxon）'}                 = qq{4_2_2_9}  if -d qq{$out_dir/Community/Wilcoxon};
			$hash{'Kruskal_Wallis差异分析（Kruskal_Wallis）'}     = qq{4_2_2_10}  if -d qq{$out_dir/Community/Kruskal_Wallis};
			$hash{'ANOVA方差分析（ANOVA）'}                       = qq{4_2_2_11} if -d qq{$out_dir/Community/ANOVA};
			$hash{'基因拷贝数矫正分析（Gene_Copies_Correction）'} = qq{4_2_2_12} if -d qq{$out_dir/Community/Gene_Copies_Correction};

			$hash{'Alpha多样性分析（AlphaDiversity)'}             = qq{4_2_3}  if -d qq{$out_dir/AlphaDiversity};
			$hash{'多样性指数统计（DiversityIndex）'}             = qq{4_2_3_1}  if -d qq{$out_dir/AlphaDiversity/DiversityIndex};
			$hash{'多样性指数差异分析（DiversityIndex）'}         = qq{4_2_3_2}  if -d qq{$out_dir/AlphaDiversity/DiversityIndex/Pairwise_groups};
			$hash{'稀释性曲线（Rarefaction）'}  			      = qq{4_2_3_3}  if -d qq{$out_dir/AlphaDiversity/Rarefaction};
			$hash{'Rank-Abundance曲线（RankAbundance）'}          = qq{4_2_3_4}  if -d qq{$out_dir/AlphaDiversity/RankAbundance};
			$hash{'物种累积曲线（Specaccum）'}                    = qq{4_2_3_5}  if -d qq{$out_dir/AlphaDiversity/Specaccum};

			$hash{'Beta多样性分析（BetaDiversity）'}              = qq{4_2_4}  if -d qq{$out_dir/BetaDiversity};
			$hash{'Venn图（Venn）'}                               = qq{4_2_4_1}  if -d qq{$out_dir/BetaDiversity/Venn};
			$hash{'样本聚类树图（SamplesTree）'}                  = qq{4_2_4_2}  if -d qq{$out_dir/BetaDiversity/SamplesTree};
			$hash{'PCA分析（PCA）'}      			              = qq{4_2_4_3}  if -d qq{$out_dir/BetaDiversity/PCA};
			$hash{'PCoA分析（PCoA）'}     			              = qq{4_2_4_4}  if -d qq{$out_dir/BetaDiversity/PCoA};
			$hash{'NMDS分析（NMDS）'}     			              = qq{4_2_4_5}  if -d qq{$out_dir/BetaDiversity/NMDS};
			$hash{'PLS-DA分析（PLS_DA）'}     			          = qq{4_2_4_6}  if -d qq{$out_dir/BetaDiversity/PLS_DA};
			$hash{'ADONIS分析（ADONIS）'}                         = qq{4_2_4_7}  if -d qq{$out_dir/BetaDiversity/ADONIS};
			$hash{'Lefse分析（Lefse）'}                           = qq{4_2_4_8}  if -d qq{$out_dir/BetaDiversity/Lefse};
			$hash{'群落构建机制分析（BetaNTI）'}                  = qq{4_2_4_9}  if -d qq{$out_dir/BetaDiversity/BetaNTI};

			$hash{'功能预测分析（PICRUSt）'}                      = qq{4_2_5}  if -d qq{$out_dir/PICRUSt};
			$hash{'功能预测'}                                     = qq{4_2_5_1}  if -d qq{$out_dir/PICRUSt};
			$hash{'单样本COG条形图（cog_barplot）'}               = qq{4_2_5_2}  if -d qq{$out_dir/PICRUSt/cog_barplot};
			$hash{'T-test分析（pair_group_different_analysis）'}  = qq{4_2_5_3}  if -d qq{$out_dir/PICRUSt/pair_group_different_analysis};
			$hash{'ANOVA分析（all_group_different_analysis）'}    = qq{4_2_5_4}  if -d qq{$out_dir/PICRUSt/all_group_different_analysis};

			$hash{'环境因子分析（Environmentfactor）'}            = qq{4_2_6}  if -d qq{$out_dir/EnvironmentFactors};
			$hash{'Mantel test分析（MantelTest）'}                = qq{4_2_6_1}  if -d qq{$out_dir/EnvironmentFactors/MantelTest};
			$hash{'Bioenv分析（BioEnv）'}                         = qq{4_2_6_2}  if -d qq{$out_dir/EnvironmentFactors/BioEnv};
			$hash{'物种与环境因子相关性分析（Correlation）'}      = qq{4_2_6_3}  if -d qq{$out_dir/EnvironmentFactors/Correlation};
			$hash{'RDA/CCA分析（RDA_CCA）'}                       = qq{4_2_6_4}  if -d qq{$out_dir/EnvironmentFactors/RDA_CCA};

		}elsif($type[$i] eq "Relative_Quantitation"){

			my $out_dir = "$out_dir1/$group_id/Relative_Quantitation";

			$hash{'相对定量数据分析（Relative_Quantitation）'}    = qq{4_3}  if -d qq{$out_dir};
			$hash{'OTU分类分析'}                                  = qq{4_3_1}  if -d qq{$out_dir/OTU};
			$hash{'OTU聚类分析'}                                  = qq{4_3_1_1}  if -d qq{$out_dir/OTU};
			$hash{'分类学分析'}                                   = qq{4_3_1_2}  if -d qq{$out_dir/OTU};
			$hash{'物种注释统计'}                                 = qq{4_3_1_3}  if -d qq{$out_dir/OTU};
			$hash{'进化树分析'}                                   = qq{4_3_1_4}  if -d qq{$out_dir/OTU};

			$hash{'群落组成分析（Community）'}                    = qq{4_3_2}  if -d qq{$out_dir/Community};
			$hash{'单样本物种组成饼图（Pieplot）'}                = qq{4_3_2_1}  if -d qq{$out_dir/Community/Pieplot};
			$hash{'单样本多级物种组成图（KronaPlot）'}            = qq{4_3_2_2}  if -d qq{$out_dir/Community/KronaPlot};
			$hash{'分类学系统组成树（TaxonTree）'}                = qq{4_3_2_3}  if -d qq{$out_dir/Community/TaxonTree};
			$hash{'物种总丰度柱状图（Community_Structure）'}      = qq{4_3_2_4}  if -d qq{$out_dir/Community/Community_Structure};
			$hash{'丰度分布Bubble图（Bubble）'}                   = qq{4_3_2_5}  if -d qq{$out_dir/Community/Bubble};
			$hash{'多样本物种组成柱状图（Barplot）'}              = qq{4_3_2_6}  if -d qq{$out_dir/Community/Barplot};
			$hash{'多样本物种组成聚类图（Treebar）'}              = qq{4_3_2_7}  if -d qq{$out_dir/Community/Treebar};
			$hash{'Heatmap热图（Heatmap）'}                       = qq{4_3_2_8}  if -d qq{$out_dir/Community/Heatmap};
			$hash{'Wilcoxon差异分析（Wilcoxon）'}                 = qq{4_3_2_9}  if -d qq{$out_dir/Community/Wilcoxon};
			$hash{'Kruskal_Wallis差异分析（Kruskal_Wallis）'}     = qq{4_3_2_10}  if -d qq{$out_dir/Community/Kruskal_Wallis};
			$hash{'Metastats分析（Metastats）'}                   = qq{4_3_2_11} if -d qq{$out_dir/Community/Metastats};
			$hash{'ANOVA方差分析（ANOVA）'}                       = qq{4_3_2_12} if -d qq{$out_dir/Community/ANOVA};

			$hash{'Alpha多样性分析（AlphaDiversity)'}             = qq{4_3_3}  if -d qq{$out_dir/AlphaDiversity};
			$hash{'多样性指数统计（DiversityIndex）'}             = qq{4_3_3_1}  if -d qq{$out_dir/AlphaDiversity/DiversityIndex};
			$hash{'多样性指数差异分析（DiversityIndex）'}         = qq{4_3_3_2}  if -d qq{$out_dir/AlphaDiversity/DiversityIndex/Pairwise_groups};
			$hash{'稀释性曲线（Rarefaction）'}  			      = qq{4_3_3_3}  if -d qq{$out_dir/AlphaDiversity/Rarefaction};
			$hash{'Shannon-Wiener曲线（Shannon）'}  			  = qq{4_3_3_4}  if -d qq{$out_dir/AlphaDiversity/Shannon};
			$hash{'Rank-Abundance曲线（RankAbundance）'}          = qq{4_3_3_5}  if -d qq{$out_dir/AlphaDiversity/RankAbundance};
			$hash{'物种累积曲线（Specaccum）'}                    = qq{4_3_3_6}  if -d qq{$out_dir/AlphaDiversity/Specaccum};


			$hash{'Beta多样性分析（BetaDiversity）'}              = qq{4_3_4}  if -d qq{$out_dir/BetaDiversity};
			$hash{'Venn图（Venn）'}                               = qq{4_3_4_1}  if -d qq{$out_dir/BetaDiversity/Venn};
			$hash{'样本聚类树图（SamplesTree）'}                  = qq{4_3_4_2}  if -d qq{$out_dir/BetaDiversity/SamplesTree};
			$hash{'PCA分析（PCA）'}      			              = qq{4_3_4_3}  if -d qq{$out_dir/BetaDiversity/PCA};
			$hash{'PCoA分析（PCoA）'}     			              = qq{4_3_4_4}  if -d qq{$out_dir/BetaDiversity/PCoA};
			$hash{'NMDS分析（NMDS）'}     			              = qq{4_3_4_5}  if -d qq{$out_dir/BetaDiversity/NMDS};
			$hash{'PLS-DA分析（PLS_DA）'}     			          = qq{4_3_4_6}  if -d qq{$out_dir/BetaDiversity/PLS_DA};
			$hash{'ADONIS分析（ADONIS）'}                         = qq{4_3_4_7}  if -d qq{$out_dir/BetaDiversity/ADONIS};
			$hash{'Lefse分析（Lefse）'}                           = qq{4_3_4_8}  if -d qq{$out_dir/BetaDiversity/Lefse};
			$hash{'群落构建机制分析（BetaNTI）'}                  = qq{4_3_4_9}  if -d qq{$out_dir/BetaDiversity/BetaNTI};

			$hash{'功能预测分析（PICRUSt）'}                      = qq{4_3_5}  if -d qq{$out_dir/PICRUSt};
			$hash{'功能预测'}                                     = qq{4_3_5_1}  if -d qq{$out_dir/PICRUSt};
			$hash{'单样本COG条形图（cog_barplot）'}               = qq{4_3_5_2}  if -d qq{$out_dir/PICRUSt/cog_barplot};
			$hash{'T-test分析（pair_group_different_analysis）'}  = qq{4_3_5_3}  if -d qq{$out_dir/PICRUSt/pair_group_different_analysis};
			$hash{'ANOVA分析（all_group_different_analysis）'}    = qq{4_3_5_4}  if -d qq{$out_dir/PICRUSt/all_group_different_analysis};

			$hash{'环境因子分析（Environmentfactor）'}            = qq{4_3_6}  if -d qq{$out_dir/EnvironmentFactors};
			$hash{'Mantel test分析（MantelTest）'}                = qq{4_3_6_1}  if -d qq{$out_dir/EnvironmentFactors/MantelTest};
			$hash{'Bioenv分析（BioEnv）'}                         = qq{4_3_6_2}  if -d qq{$out_dir/EnvironmentFactors/BioEnv};
			$hash{'物种与环境因子相关性分析（Correlation）'}      = qq{4_3_6_3}  if -d qq{$out_dir/EnvironmentFactors/Correlation};
			$hash{'RDA/CCA分析（RDA_CCA）'}                       = qq{4_3_6_4}  if -d qq{$out_dir/EnvironmentFactors/RDA_CCA};

		}elsif($type[$i] eq "Combine_Analysis"){

			my $out_dir = "$out_dir1/$group_id/Combine_Analysis";

			$hash{'联合分析（Combine_Analysis）'}                 = qq{4_4}  if -d qq{$out_dir};	
			$hash{'联合分析简介'}                                 = qq{4_4_1}  if -d qq{$out_dir};	

			$hash{'群落组成对比分析（Community）'}                = qq{4_4_2}  if -d qq{$out_dir/Community};
			$hash{'多样本物种组成柱状图对比分析（Barplot）'}      = qq{4_4_2_1}  if -d qq{$out_dir/Community/Barplot};
			$hash{'Wilcoxon结果比对分析（Wilcoxon）'}             = qq{4_4_2_2}  if -d qq{$out_dir/Community/Wilcoxon};
			$hash{'Kruskal_Wallis结果对比分析（Kruskal_Wallis）'} = qq{4_4_2_3}  if -d qq{$out_dir/Community/Kruskal_Wallis};
			$hash{'ANOVA结果对比分析（ANOVA）'}                   = qq{4_4_2_4}  if -d qq{$out_dir/Community/ANOVA};

			$hash{'Alpha多样性对比分析（AlphaDiversity）'}        = qq{4_4_3}  if -d qq{$out_dir/AlphaDiversity};
			$hash{'多样性指数对比分析（DiversityIndex）'}         = qq{4_4_3_1}  if -d qq{$out_dir/AlphaDiversity/DiversityIndex};

			$hash{'Beta多样性对比分析（BetaDiversity）'}          = qq{4_4_4}  if -d qq{$out_dir/BetaDiversity};
			$hash{'ADONIS结果对比分析（ADONIS）'}                 = qq{4_4_4_1}  if -d qq{$out_dir/BetaDiversity/ADONIS};
			$hash{'Lefse结果对比分析（Lefse）'}                   = qq{4_4_4_2}  if -d qq{$out_dir/BetaDiversity/Lefse};

			$hash{'功能预测对比分析（PICRUSt）'}                  = qq{4_4_5}  if -d qq{$out_dir/PICRUSt};
			$hash{'T-test结果对比分析（pair_group_different_analysis）'} = qq{4_4_5_1}  if -d qq{$out_dir/PICRUSt/pair_group_different_analysis};
			$hash{'ANOVA结构对比分析（all_group_different_analysis）'}   = qq{4_4_5_2}  if -d qq{$out_dir/PICRUSt/all_group_different_analysis};			

		}
		
	}

	return \%hash;
}


sub change_catalog_stepid {

	my ($catalog_tmp, $catalog) = @_;

	my @id = split / /, `cat $catalog_tmp | cut -f 2 |xargs`;
	my @step = split / /, `cat $catalog_tmp | cut -f 1 |xargs`;

	for (my $i = 0; $i < @id; $i++) {
		
		if ($id[$i] =~ /^4/) {
			my $f = $id[$i-1];
			my $n = $id[$i+1];

			if (($f =~ /^\d+_\d+_\d+_\d+$/ or $f =~ /^\d+_\d+_\d+$/) and $n =~ /^\d+_\d+_\d+$/) {
				
				$step[$i] = "S10" if ( $step[$i] ne "S9"); 

			}elsif($f =~ /^\d+_\d+_\d+_\d+$/ and $n =~ /^\d+_\d+$/){
				$step[$i] = "S13";

			}elsif($n =~ /^5/){

				$step[$i] = "S14";
			}

		}
	}


	open CATA, "$catalog_tmp";
	open CATALOG,">$catalog";
	my $cnt = 0;
	while (<CATA>){
		chomp;
		chomp $step[$cnt];
		my ($id, $tmp) = split /\t/, $_, 2;
		print CATALOG qq{$step[$cnt]\t$tmp\n};	
		$cnt++;

	}

	close CATA;
	close CATALOG;

}

1
