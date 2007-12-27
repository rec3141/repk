#!/opt/local/bin/perl
###
### Switched orientation of enzmatrix in v.0.3.0 so need to recorrect it here
###

 use strict;
# use lib '/usr/lib/perl5/site_perl/5.8.0';
# use lib '/da41/d61/rec3141/public_html/lib/';
# use lib '~/public_html/lib';
# use lib '/opt/local/lib/perl5/vendor_perl/5.8.8/darwin-2level/';

 use GD::Simple;

my ($matrix_ref, $zymes_ref, $names_ref) = readFile(@ARGV);
my @matrix = @$matrix_ref;
my @names = @$names_ref;
my @zymes = @$zymes_ref;

# print join(',',@names);
# print join(',',@zymes);

#     for my $i ( 0 .. $#matrix ) {
#         print "\trow $i ($names[$i]) is [ @{$matrix[$i]} ],\n";
#     }

my @species_names = @names;
my $rows = @names;
my $cols = @zymes;

#  my @color_names = GD::Simple->color_names;
#  my $cols = int(sqrt(@color_names));
#  my $rows = int(@color_names/$cols)+1;



 my $cell_width    = 20;
 my $cell_height   = 20;
 my $spacer_height = 3;
 my $spacer_width = 3;
 my $legend_width = 250;
 my $title_height = 100;
 my $width       = $cols * $cell_width + 2*$legend_width;
 my $height      = $rows * $cell_height + $title_height;

 my $img = GD::Simple->new($width,$height);
 $img->font(gdSmallFont);

my @topleft;
my @botright;# = (0, $origin[1]); #+cell_height
my @textsite = (0, $title_height+$cell_height);
my $species;
my $enzyme;

 for (my $r=0; $r<$rows; $r++) {
     $species = $species_names[$r] or next;
     @textsite = (0, ($r+1)*$cell_height + $title_height - $spacer_height);
     $img->moveTo(@textsite);
     $img->fgcolor('black');
     $img->angle(0);
     $img->string($species);

   for (my $c=0; $c<$cols; $c++) {
     @topleft  = ($c*$cell_width+$legend_width,$r*$cell_height+$title_height);
     @botright = ($topleft[0]+$cell_width-$spacer_width,$topleft[1]+$cell_height-$spacer_height);
     $img->fgcolor('blue');
     if ($matrix[$r][$c] == 1) {$img->bgcolor('blue')}
     else {$img->bgcolor('white')}
# $img->moveTo(@topleft);
# $img->string("$c-$r");
     $img->rectangle(@topleft,@botright);

#print enzyme names during each column
     $enzyme = $zymes[$c] or next;
     $img->moveTo($topleft[0]+0.8*$cell_width,$title_height);
     $img->fgcolor('black');
     $img->angle(-90);
     $img->string($enzyme);
   }
 }

 print $img->png;


sub readFile{
	my $inputFile = shift; #'./file.csv'; 
	my @matrix;
	my @names;
	my @zymes;
	my $count = 1;
	open  (FH, $inputFile) or die "Couldn't open $inputFile: $!";
	while (my $line =<FH>) {
# 		$line =~ s/\r//;
		chomp $line;
		my @array = split (/\s/, $line);   #changed from \t to \s
# 		print join(',',@array) . "\n";
		if ($count == 1) {
			$count++;
			shift @array;
			@names = @array;
			
		}

		else{
			my $enzyme = shift(@array);
			push @zymes, $enzyme;
			push @matrix, [ @array ];
		}
	}	

	close FH or die "Couldn't close $inputFile: $!";
	return(\@matrix,\@names, \@zymes);
}
