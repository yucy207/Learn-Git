#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($node, $edge, $out_dir,  $help);

GetOptions(
    'node|n=s'    => \$node,
    'edge|e=s'    => \$edge,
    'out_dir|o=s' => \$out_dir,
    'help|h!'     => \$help
);

if ($help or not $node or not $edge or not $out_dir) {
    usage();
    exit();
}

sub usage
{
my $help =<<EOF;
Usage: perl cytoscape_js.pl -node nodes.txt -edge edges.txt -out_dir out_dir
    -node     -n  nodes txt  file [required]
    -edge     -e  edges txt  file [required]
    -out_dir  -o  output dir [required]
    -help     -h  print help message
EOF

print $help;

}

my %nodes_attr = ();
my %edges_attr = ();


download_js();
parse_nodes();
parse_edges();
write_html();
html_2_png();

sub download_js
{

    system qq{mkdir -p $out_dir/js} if not -d qq{$out_dir/js};
    system qq{wget -q -O $out_dir/js/jquery.min.js     "http://cdn.bootcss.com/jquery/1.11.2/jquery.min.js"};
    system qq{wget -q -O $out_dir/js/cytoscape.min.js  "http://cdn.bootcss.com/cytoscape/2.3.16/cytoscape.min.js"};
}


sub parse_nodes
{

    my $cnt  = 0;
    open TXT, $node or die "Can't open $node!\n";
    while (<TXT>) {
        chomp;
        $cnt++;
        my @arr = split /\t/;

        $nodes_attr{$arr[0]} = [$cnt, $arr[1], $arr[2]];
    }
    close TXT;
}

sub parse_edges
{
    open TXT, $edge or die "Can't open $edge!\n";
    while (<TXT>) {
        chomp;
        $edges_attr{$_} = "";
    }
    close TXT; 
}

sub write_html
{

open SAVE, qq{>$out_dir/cytoscape.html};
my @datas = ();

foreach my $x (keys %nodes_attr) {
    my $id   = $nodes_attr{$x}->[0];
    my $type = $nodes_attr{$x}->[1];
    my $regulation = $nodes_attr{$x}->[2];

    my $line = qq{               {data : {id : '$id', name : '$x', type : '$type', regulation : '$regulation'}}};
    push @datas, $line;
}

my $nodes_info = join ",\n", @datas;




my @edge_datas = ();
foreach my $x (keys %edges_attr) {
    my ($source, $target) = split /\t/, $x;
    my $source_id   = $nodes_attr{$source}->[0];
    my $target_id   = $nodes_attr{$target}->[0];


    my $line = qq{              {data : {source : '$source_id', target : '$target_id'}}};
    push @edge_datas, $line;
}
my $edges_info = join ",\n", @edge_datas;




my $help =<<EOF;
<!DOCTYPE html>
<html>
<head>
    <title>Learning Cytoscape.js</title>
    <style type="text/css">
        /* cytoscape graph */
        #cy {
            height: 1200px;
            width: 1800px;
            background-color: #f9f9f9;
        }
    </style>
    <script src="js/jquery.min.js"></script>
    <script src="js/cytoscape.min.js"></script>
    <script>
        \$(function(){
            cytoscape({
              container: document.getElementById('cy'),
              style: [
                { selector: 'node[type = "mRNA"]', 
                  css: { 'content': 'data(name)', 'shape': 'triangle'}
                },
                { selector: 'node[type = "miRNA"]', 
                  css: { 'content': 'data(name)', 'shape': 'rectangle'}
                },
                { selector: 'node[type = "circRNA"]', 
                  css: { 'content': 'data(name)', 'shape': 'ellipse'}
                },
                { selector: 'node[type = "lncRNA"]', 
                  css: { 'content': 'data(name)', 'shape': 'ellipse'}
                },
                { selector: 'node[regulation = "Up"]', 
                  css: {'background-color': 'red'}
                },
                { selector: 'node[regulation = "Down"]', 
                  css: {'background-color': 'green'}
                }                                        
              ],
              elements: {
                nodes:[\n$nodes_info\n],
                edges:[\n$edges_info\n]
              },
              layout: { name: 'cose'} 
            });
        });
    </script>
</head>
<body>
    <div id="cy"></div>
</body>
</html>  
EOF

print SAVE $help;
close SAVE;
}

sub html_2_png
{
my $txt =<<EOF;
var page = require('webpage').create();
page.open('$out_dir/cytoscape.html', function(status) {
    setTimeout(function(){
        page.render('$out_dir/network.png');
        phantom.exit();
    }, 10000);
});
EOF

open SAVE, qq{>$out_dir/js/html2png.js};
print SAVE qq{$txt\n};
close SAVE;

system qq{/home/genesky/software/phantomjs/2.1.1/bin/phantomjs $out_dir/js/html2png.js\n};
}