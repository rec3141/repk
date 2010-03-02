#!/opt/local/bin/perl
# -w
# this program downloads and parses the lastest version of the REBASE database
# for use in REPK (http://staff.washington.edu/rec3141/repk)
# this is version 2

use strict;

# use LWP::Simple;
# use Data::Dumper;

# data source from REBASE:
# 8. Type II (tabbed output)

#my $url = 'ftp://ftp.neb.com/rebase/link_itype2';
#my $rebase = get($url);
#print "blah $rebase";

my $repk_location="/Library/Webserver/Documents/cgi/repk/";
my $wget_location="/sw/bin/wget";
`rm link_itype2`;
`$wget_location --user-agent="Mozilla" http://rebase.neb.com/rebase/link_itype2`;

my %Enzymes;
my $version;
my ($enzyme, $prototype, $recsite, $junk, $comsource);
open(FILEIN, "< link_itype2");

LINE: while (<FILEIN>) {

if (m/version/) {
 my @versionline = split('\s');
 $version = $versionline[2];
}

if (m/\t/) {
 ($enzyme, $prototype, $recsite, $junk, $comsource) = split("\t");
 chomp ($enzyme, $prototype, $recsite, $comsource);
}
else {next LINE};

if (length($comsource) > 0) {}
else {next LINE}

# next LINE unless (length($comsource) == 111);

   if (length($Enzymes{$recsite}{'prototype'}) > 0) {
    $Enzymes{$recsite}{'isos'} .= "$enzyme\t";
    }
   else {
    $Enzymes{$recsite}{'prototype'} = "$enzyme";
   }
}


close FILEIN or die "complete failure";



#now need to sort out all those that have () in their names, although will need to deal with these eventually.
my %NewEnzymes;
foreach my $enzkey (keys %Enzymes) {
my $temptype = $Enzymes{$enzkey}{'prototype'};
$NewEnzymes{$temptype}{'isos'} = $Enzymes{$enzkey}{'isos'};
$NewEnzymes{$temptype}{'recsite'} = $enzkey;


}

my %WeirdEnzymes;
foreach my $enzkey2 (keys %NewEnzymes) {
 if ($NewEnzymes{$enzkey2}{'recsite'} =~ m/\(/) {
  $WeirdEnzymes{$enzkey2} = $NewEnzymes{$enzkey2};
  delete $NewEnzymes{$enzkey2};
 }
}

KEY: foreach my $enzkey7 (keys %WeirdEnzymes) {
# next if $WeirdEnzymes{$enzkey7}{'recsite'} =~ m/^\(/;
if ($WeirdEnzymes{$enzkey7}{'recsite'} =~ m/^\(/) {
 delete $WeirdEnzymes{$enzkey7};
next KEY}

$WeirdEnzymes{$enzkey7}{'recsite'} =~ m/^([a-zA-Z]*)\(([0-9-]{1,2})\/([0-9-]{1,2})\)/;
if ($2 > 0) {$WeirdEnzymes{$enzkey7}{'recsite'} = "$1" . ("N" x $2) . "^";}
else {$WeirdEnzymes{$enzkey7}{'recsite'} = substr($1,0,$2) . "^" . substr($1,$2)}

    ##### reverse complements each sequence to search for opposite strand cutsites
#     my $new1 = reverse( $1 );
    my $new1 = $1;
    $new1 =~ tr/ACGTUMRWSYKVHDBXN/TGCAAKYWSRMBDHVXN/;


if ($3 > 0) {$WeirdEnzymes{$enzkey7}{'altsite'} = "^" . reverse("$new1" . ("N" x $3));}
else {$WeirdEnzymes{$enzkey7}{'altsite'} = reverse(substr($new1,0,$3) . "^" . substr($new1,$3))}
}


open(OUT, "> $repk_location"."enzymes_type2p.txt");
foreach my $enzkey3 (sort {$a cmp $b} keys %NewEnzymes) {
print OUT "$enzkey3\t$NewEnzymes{$enzkey3}{'recsite'}\tISO\t$NewEnzymes{$enzkey3}{'isos'}\n";
}

#foreach my $enzkey4 (sort {$a cmp $b} keys %WeirdEnzymes) {
#print OUT "$enzkey4\t$WeirdEnzymes{$enzkey4}{'recsite'}\tISO\t$WeirdEnzymes{$enzkey4}{'isos'}\n";
#print OUT "$enzkey4+\t$WeirdEnzymes{$enzkey4}{'altsite'}\tISO\t$WeirdEnzymes{$enzkey4}{'isos'}\n";
#}
print OUT "NOCUT\tZZZ^ZZZ\tISO\tFakeEnzyme";
close OUT or die "complete failure";

open(OUT, "> $repk_location"."enzymes_type2a.txt");
foreach my $enzkey4 (sort {$a cmp $b} keys %WeirdEnzymes) {
print OUT "$enzkey4\t$WeirdEnzymes{$enzkey4}{'recsite'}\tISO\t$WeirdEnzymes{$enzkey4}{'isos'}\n";
print OUT "$enzkey4+\t$WeirdEnzymes{$enzkey4}{'altsite'}\tISO\t$WeirdEnzymes{$enzkey4}{'isos'}\n";

}
close OUT or die "complete failure";

open(OUT, "> $repk_location"."rebase_version.txt");
print OUT "$version";
close OUT or die "complete failure";




# 
# # print Dumper %NewEnzymes;
# # 
# open(OUT, "> enzymes_type2_testing2.txt");
# foreach my $enzkey6 (sort {
# 	($a =~ m/([a-zA-Z]*)\^([a-zA-Z]*)/); my $aa = $1 . $2;
# 	($b =~ m/([a-zA-Z]*)\^([a-zA-Z]*)/); my $bb = $1 . $2;
# 	$aa cmp $bb} 
# 	keys %Enzymes) {
# 
# print OUT "$enzkey6\t$Enzymes{$enzkey6}{'prototype'}\tISO\t$Enzymes{$enzkey6}{'isos'}\n";
# }
# print OUT "NOCUT\tZZZ^ZZZ\tISO\tFakeEnzyme";
# close OUT or die "complete failure";
# 
# `sort -d enzymes_type2_testing2.txt > enzymes_type2_testing2a.txt`;
